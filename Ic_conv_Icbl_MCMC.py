from __future__ import print_function

# -*- coding: utf-8 -*-
##############################################################################
# Version 1.0, July 2016
'''
Fits a SN template spectrum to a SN observed spectrum 
to measure spectral features
blue-shift and broadening


 Arguments: 
 flattened Ic-bl spectrum (use snidflat.pro) in .save IDL format
 the corresponding uncertainty array
 SNe Ic template
 Output:
 marginalized distribution of model parameters, including but not
 limited to absorption velocity, see Modjaz et al. (2016) for details
 Note: initial values and prior of model parameters, and region to find
 initial template fitting region can be changed as needed
'''
import sys
import time
import corner
import emcee
import pickle as pkl
import numpy as np
import pylab as pl
from matplotlib import  rcParams
from matplotlib.ticker import MultipleLocator
from scipy.ndimage import filters
from scipy.signal import gaussian
from scipy.interpolate import interp1d
from scipy.io.idl import readsav


# initial parameter values
# v/1000 in km/s, sigma/10 in angstrom, y-amplitude, wave-range in angstrom
P0Fe = np.array([11.0, 1.0, 1.0, 1.0])
# prior for v/1000, sigma/10, y-amplitude
PriorFe = np.array([0, 30, 0, 5, 0, 3])
# region to find initial template fit region
X0Fe = np.array([4200, 4800, 5000, 5600])


def readdata(spec, template):
    ''' readdata: read in flattened Ic-bl spectrum and the corresponding
    uncertainty array,  SNe Ic template '''

    # read in Ic template
    try:
        s = readsav(template)
    except  Exception:
        print("readsav filed. Must pass 2 .sav files as input")
        return [-1]*6
    wlog_input = s.wlog[np.where((s.wlog > 4400) & (s.wlog < 9000))]
    fmean_input = s.fmean[np.where((s.wlog > 4400) & (s.wlog < 9000))]

    # read in Icbl spectrum and uncertainty array
    s2 = readsav(spec)
    x_flat = s2.wavelog_input[0:1024]
    y_flat_sm = s2.flatflux_input_sm
    y_flat = s2.flatflux_input
    y_flat_err = s2.flatflux_err_input

    return wlog_input, fmean_input, x_flat, y_flat_sm, y_flat, y_flat_err


def fittemplate(p, fmean_input, wlog_input, lx, ly, ly_err, x_flat, y_flat,
                ax=None):
    ''' fittemplate: fit SN Ic-bl spectrum to a convolved and blueshifted SN Ic template'''

    # assign values to parameters of absorption velocity, width of absorption
    # feature, amplitude, and wavelength-range
    v, sig, amplitude, w_range = -p[0] * 1000, p[1] * 10, p[2], p[3]

    # template fit region
    w_lower = lx[0] + w_range
    w_upper = lx[-1] - w_range
    inds = (lx > w_lower) * (lx < w_upper)
    ly_new = ly[inds]
    ly_err_new = ly_err[inds]
    lx_new = lx[inds]

    # Gaussian convolution
    b = gaussian(300, sig)
    thisy = filters.convolve1d(amplitude * fmean_input, b / b.sum())
    # blue shifted
    beta = v / 299792.458
    doppler = np.sqrt((1 + beta) / (1 - beta))
    f2 = interp1d(wlog_input * doppler, thisy, bounds_error=False,
                  fill_value=0)(lx_new)
    chisq = np.sum((ly_new - f2) ** 2 / ly_err_new ** 2) / (len(ly_new) - len(p))

    if not ax is None:
        minorLocatory = MultipleLocator(0.02)
        minorLocatorx = MultipleLocator(100)
        ax.xaxis.set_minor_locator(minorLocatorx)
        ax.yaxis.set_minor_locator(minorLocatory)

        ax.plot(x_flat, y_flat, 'k', alpha=0.5, label="Input SN Ic-bl spectrum")
        ax.plot(lx_new, ly_new, 'k', linewidth=3, label="Fitting region")
        ax.plot(wlog_input * doppler, amplitude * fmean_input, 'r', linewidth=2,
                alpha=0.5, label="Blueshifted SN Ic template")
        ax.plot(lx_new, f2, 'r', linewidth=4, label="Fitted SN Ic template")
        ax.text(4500, -0.4, r"$v$=%.0f km s$^{-1}$   $\sigma$=%.0f km s$^{-1}$" % 
                (-v, sig * 400), fontsize=25)
        ax.text(4500, -0.5, r"$a$=%.1f   $\Delta$$w$=%.0f $\AA$" %
                (amplitude, w_range), fontsize=25)
        ax.text(5500, 0.3, r"$\chi^2_r$=%.1f" % (chisq), fontsize=25)
        ax.set_xlabel("Rest Wavelength ($\AA$)", fontsize=25)
        ax.set_ylabel("Relative Flux", fontsize=25)
        ax.legend(fontsize=25)

    return chisq


def logprior(p):
    v = p[0]
    s = p[1]
    amplitude = p[2]
    w_range = p[3]
    if s > PriorFe[2] and s < PriorFe[3] and v > PriorFe[0] and \
       v < PriorFe[1] and amplitude > PriorFe[4] and amplitude < PriorFe[5]:
        return np.log(np.exp(-(w_range) ** 2 / (2 * 33 ** 2)))
    return -np.inf


def logl(p, x, y, s, fmean_input, wlog_input, x_flat, y_flat):
    # log likelihood
    return -np.log(s) - 0.5 * (fittemplate(p, fmean_input, wlog_input,
                                           x, y, s, x_flat, y_flat))


def logp(p, x, y, s, fmean_input, wlog_input, x_flat, y_flat):
    # full log probability 
    lgl = logl(p, x, y, s, fmean_input, wlog_input, x_flat, y_flat)
    return np.sum(logprior(p) + lgl)


def runMCMC(wlog_input, fmean_input, x_flat, y_flat_sm, y_flat, y_flat_err,
            spec, x0=X0Fe, p00=P0Fe, prior=PriorFe, posterior_save=False,
            plot_save=False, file_save=False, plotChain=False):
    '''runMCMC: sample probability distribution using package emcee,
    and get marginalized distribution of model parameters  '''

    ndim, nwalkers = 4, 6 * 2

    # template fit region
    Fe_lower_inds = (x_flat > x0[0]) * (x_flat < x0[1])
    Fe_lower_x_flat = x_flat[Fe_lower_inds]
    Fe_lower_y_flat_sm = y_flat_sm[Fe_lower_inds]
    Fe_lower = np.min(Fe_lower_x_flat[np.where(Fe_lower_y_flat_sm ==
                      np.max(Fe_lower_y_flat_sm))]) + 100
    Fe_upper_inds = (x_flat > x0[2]) * (x_flat < x0[3])
    Fe_upper_x_flat = x_flat[Fe_upper_inds]
    Fe_upper_y_flat_sm = y_flat_sm[Fe_upper_inds]
    Fe_upper = np.max(Fe_upper_x_flat[np.where(Fe_upper_y_flat_sm ==
                      np.max(Fe_upper_y_flat_sm))]) - 100

    Fex = x_flat[np.where((x_flat > Fe_lower) & (x_flat < Fe_upper))]
    Fey = y_flat[np.where((x_flat > Fe_lower) & (x_flat < Fe_upper))]
    Fes = y_flat_err[np.where((x_flat > Fe_lower) & (x_flat < Fe_upper))]

    best_pos = []
    p0 = [p00 + 1e-6 * np.random.randn(ndim) for i in range(nwalkers)]
    samplerFe = emcee.EnsembleSampler(nwalkers, ndim, logp,
                                      args=(Fex, Fey, Fes, fmean_input,
                                            wlog_input, x_flat, y_flat))
    # run MCMC for 30 steps starting from the tiny ball defined above
    samplerFe.run_mcmc(p0, 30)
    best_pos.append(samplerFe.flatchain[samplerFe.flatlnprobability.argmax()])
    pos = emcee.utils.sample_ball(best_pos[-1], best_pos[-1] / 100.,
                                  size=nwalkers)
    samplerFe.reset()
    samplerFe.run_mcmc(pos, 1000)
    best_pos.append(samplerFe.flatchain[samplerFe.flatlnprobability.argmax()])
    value_50 = [np.percentile(samplerFe.chain[:
        , :
        , 0], [50])[0],
                np.percentile(samplerFe.chain[:
        , :
        , 1], [50])[0],
                np.percentile(samplerFe.chain[:
        , :
        , 2], [50])[0],
                np.percentile(samplerFe.chain[:
        , :
        , 3], [50])[0]]

    # save marginalized distribution of model parameters
    if posterior_save:
        with open(posterior_save, 'w') as f:
            pkl.dump(samplerFe, f)

    # save template fit plot, corner plot (and chain plot if requested)
    if plot_save:

        # save template fit plot
        fig, ax = pl.subplots(figsize=(15, 15))
        
        fittemplate(value_50, fmean_input, wlog_input, Fex, Fey, Fes, x_flat,
                    y_flat, ax = ax)
        
        pl.savefig(plot_save)
        pl.close(fig)

        # save corner plot
        value_50 = [np.percentile(samplerFe.chain[:
            , :
            , 0], [50])[0],
                    np.percentile(samplerFe.chain[:
            , :
            , 1], [50])[0] * 4,
                    np.percentile(samplerFe.chain[:
            , :
            , 2], [50])[0],
                    np.percentile(samplerFe.chain[:
            , :
            , 3], [50])[0]]

        samplerFe.flatchain[:
            , 1] = samplerFe.flatchain[:
            , 1] * 4

        rcParams['xtick.labelsize'] = 22.
        rcParams['ytick.labelsize'] = 22.
        fig_corner = corner.corner(samplerFe.flatchain, truths=value_50,
                      quantiles=[0.16, 0.5, 0.84],
                      labels=[r"$v$ [10$^3$ km s$^{-1}$]",
                              "$\sigma$ [$10^3$ km s$^{-1}$]", "$a$",
                              "$\Delta$$w$ [$\AA$]"])
        fig_corner.savefig(plot_save.replace('Fe.pdf','FeFit.pdf'))
        pl.close(fig_corner)

        # save chain plot
        if plotChain:
            y_label = [r"$v/1000$ [km s$^{-1}$]", "$\sigma/1000$ [km s$^{-1}$ ]",
                       "amplitude", "wave-range [$\AA$]"]
            figChain = pl.figure(figsize=(15, 4*ndim))
            for i in range(ndim):
                ax = figChain.add_subplot(ndim,1,i+1)
                ax.plot(range(1000), samplerFe.chain[:
                                                     , :
                                                     , i].T)
                if i == ndim - 1:
                    ax.set_xlabel("steps", fontsize=30)
                ax.set_ylabel(y_label[i], fontsize=30)
            figChain.savefig(plot_save.replace('Fe.pdf','FeChain.pdf'))
        pl.close('all')


    # save initial template fit region, mean acceptance fraction, initial
    # values for parameters, and 16th, 50th, 84th percentiles of marginalized
    # distribution of model parameters
    if file_save:
        f = open(file_save, 'w')
        f.write('region to find initial template fit region:' + str(x0) + '\n')
        f.write("Mean acceptance fraction: {0:.3f} \n"
                .format(np.mean(samplerFe.acceptance_fraction)))
        f.write('uniform prior for v/1000 in km/s, sigma/10 in angstrom, ' +
                'amplitude: ' + str(prior) + '\n')
        f.write('v/1000 in km/s, sigma/1000 in km/s, amplitude, ' +
                'wave-range in angstrom' + '\n')
        f.write('initial guess: ' + str(p00) + '\n')
        f.write('best value: ' + str(best_pos[-1]) + '\n')
        f.write('16th, 50th, 84th percentiles \n')
        f.write(str(np.percentile(samplerFe.chain[:
            , :
            , 0], [16, 50, 84])) +
                ' for v/1000 in km/s\n')
        f.write(str(np.percentile(samplerFe.chain[:
            , :
            , 1], [16, 50, 84])) +
                ' for sigma/1000 in km/s\n')
        f.write(str(np.percentile(samplerFe.chain[:
            , :
            , 2], [16, 50, 84])) +
                ' for amplitude\n')
        f.write(str(np.percentile(samplerFe.chain[:
            , :
            , 3], [16, 50, 84])) +
                ' for wave-range in angstrom\n')
        f.close()

    print("Mean acceptance fraction: {0:.3f}"
          .format(np.mean(samplerFe.acceptance_fraction)))
    print('16th, 50th, 84th percentiles of marginalized distribution' +
          'of model parameters')
    # 16th, 50th, 84th percentiles of the velocity/1000
    print(str(np.percentile(samplerFe.chain[:
        , :
        , 0], [16, 50, 84])) + \
        ' for v/1000 in km/s')
    # 16th, 50th, 84th percentiles of the sigma/10000 in km/s
    print(str(np.percentile(samplerFe.chain[:
        , :
        , 1], [16, 50, 84])) + \
        ' for sigma/1000 in km/s')
    # 16th, 50th, 84th percentiles of the amplitude
    print(str(np.percentile(samplerFe.chain[:
        , :
        , 2], [16, 50, 84])) + \
        ' for amplitude')
    # 16th, 50th, 84th percentiles of the wave-range in angstrom
    print(str(np.percentile(samplerFe.chain[:
        , :
        , 3], [16, 50, 84])) + \
        ' for wave-range in angstrom')


def conv(spec, template):

    try:
        wlog_input, fmean_input, x_flat, y_flat_sm, y_flat, y_flat_err = \
                                                    readdata(spec, template)
    except IOError:
        print("must pass 2 .save files as input")
        return -1
    if isinstance(wlog_input, int):
        return -1
    if np.mean(y_flat[np.where(x_flat < 4500)]) != 0 and \
       np.mean(y_flat[np.where(x_flat > 5100)]) != 0:
        runMCMC(wlog_input, fmean_input, x_flat, y_flat_sm,
                y_flat, y_flat_err, spec,
                posterior_save=spec.replace('.sav','') + '-Fe.p',
                plot_save=spec.replace('.sav','') + '-Fe.pdf',
                file_save=spec.replace('.sav','') + '-Fe.dat', plotChain=True)
    else:
        print("wavelength range doesn't match")
        return -1
    return 0
if __name__ == '__main__':

    helpStrg = '''Arguments: observed spectrum and template spectrum: use as
$python Ic_conv_Icbl_MCMC.py 10qts_20100815_Lick_3-m_v1-z.flm-flat.sav  meanspecIc_0.sav'''
    
    t1 = time.time()

    if len(sys.argv) == 1:
        # use default arguments for testing
        obsSpec = '10qts_20100815_Lick_3-m_v1-z.flm-flat.sav'
        templSpec = 'meanspecIc_0.sav'

    elif len(sys.argv) == 3:
        # assume that you passed correctly the spectrum and template spectrum
        obsSpec, templSpec = sys.argv[1:]

    else:
        print(helpStrg)
        sys.exit()

    print("running convolution...")
    if conv(obsSpec, templSpec) == -1:
        print(helpStrg)
        sys.exit()

    t2 = time.time()
    print('minimization took {} seconds'.format(t2 - t1))
