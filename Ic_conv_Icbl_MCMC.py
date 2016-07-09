from __future__ import print_function

# -*- coding: utf-8 -*-
##############################################################################
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
import os
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

# output directory
outdir = "./outputs"

# spectra range allowed: we restrict the usable wavelength range
# generally the quality of IR spectra and the Fe blanketing in UV limit
# the region of the spectrum where our method should be used
Wlim = (4400, 9000)

from elementDicts import *


def readdata(spec, template):
    ''' readdata: read in flattened Ic-bl spectrum and the corresponding
    uncertainty array,  SNe Ic template '''

    print("reading inputs...")

    # read in Ic template
    try:
        try:
            phase = int(float(template))
            template = ("IcTemplates/meanspecIc_%d.sav"%phase)
            s = readsav(template)            
        except ValueError:
            try:
                phase = int(template.split("_")[-1].replace(".sav",""))
            except ValueError:
                phase = np.nan
            s = readsav(template)
    except IOError:
        
        print("\nError: Failing while reading the template", template)
        print("You must pass 2 files as input: ")
        print("- a spectrum file in .sav or .csv format")
        print("- a template spectrum file in .sav format ")
        print("  or a phase (number of days since Vmax, min=-10 max=72) if using the meantemplate distributed with this package\n")        
        return [-1]*6
    
    # reads in spectrum        
    try:
        if spec.endswith('.sav'):
            s2 = readsav(spec)
        else:
            try:
                # read csv in with pandas if installed
                import pandas as pd
                s2 = pd.read_csv(spec).to_dict(orient="list")
            except ImportError:
                # read csv with numpy and convert it to dictionary
                tmp = np.genfromtxt(spec, delimiter=',', dtype=None)
                s2 = {}
                s2['wavelog_input'] = tmp[1:,0].astype(float)                
                s2['flatflux_input'] = tmp[1:,1].astype(float) 
                s2['flatflux_err_input'] = tmp[1:,2].astype(float) 
                s2['flatflux_input_sm'] = tmp[1:,3].astype(float) 
                
    except  Exception:
        print("\nFailing while reading the input spectrum", spec)
        print("You must pass 2 files as input: ")
        print("- a spectrum file in .sav or .csv format")
        print("- a template spectrum file in .sav format ")
        print("  or a phase (integer number of days) if using the meantemplate distributed with this package\n")        
        return [-1]*6

    # we restrict the range to 4400 9000 A where most spectra are well sampled
    wlog_input = s['wlog'][(s['wlog'] > Wlim[0]) * (s['wlog'] < Wlim[1])]
    fmean_input = s.fmean[(s['wlog'] > Wlim[0]) * (s['wlog'] < Wlim[1])]

    # read in Icbl spectrum and uncertainty array

    y_flat_sm = np.array(s2['flatflux_input_sm'])
    y_flat = np.array(s2['flatflux_input'])
    y_flat_err = np.array(s2['flatflux_err_input'])
    
    # the default wavelength array created from our templates
    # has one extra value at the end
    x_flat = np.array(s2['wavelog_input'][:y_flat.size])

    return wlog_input, fmean_input, x_flat, y_flat_sm, y_flat, y_flat_err, phase


def fittemplate(p, element, fmean_input, wlog_input,
                lx, ly, ly_err, x_flat, y_flat,
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
        ax.plot(x_flat, y_flat, 'k', alpha=0.5,
                label="Input SN Ic-bl spectrum")
        ax.plot(lx_new, ly_new, 'k', linewidth=3, label="Fitting region")
        ax.plot(wlog_input * doppler, amplitude * fmean_input, 'r',
                linewidth=2, alpha=0.5, label="Blueshifted SN Ic template")
        ax.plot(lx_new, f2, 'r', linewidth=4, label="Fitted SN Ic template")
        ax.text(4500, -0.4, r"$v$=%.0f km s$^{-1}$   $\sigma$=%.0f km s$^{-1}$"
                % (-v, sig * 400), fontsize=25)
        ax.text(4500, -0.5, r"$a$=%.1f   $\Delta$$w$=%.0f $\AA$" %
                (amplitude, w_range), fontsize=25)
        ax.text(5500, 0.3, r"$\chi^2_r$=%.1f" % (chisq), fontsize=25)
        ax.set_xlabel("Rest Wavelength ($\AA$)", fontsize=25)
        ax.set_ylabel("Relative Flux", fontsize=25)
        ax.legend(fontsize=25)

    return chisq


def logprior(p, element):
    v = p[0]
    s = p[1]
    amplitude = p[2]
    w_range = p[3]
    if s > Prior[element][2] and s < Prior[element][3] and \
       v > Prior[element][0] and v < Prior[element][1] and \
       amplitude > Prior[element][4] and amplitude < Prior[element][5]:
        return np.log(np.exp(-(w_range) ** 2 / (2 * 33 ** 2)))
    return -np.inf


def logl(p, element, x, y, s, fmean_input, wlog_input, x_flat, y_flat):
    # log likelihood
    return -np.log(s) - 0.5 * (fittemplate(p, element, fmean_input, wlog_input,
                                           x, y, s, x_flat, y_flat))


def logp(p, element, x, y, s, fmean_input, wlog_input, x_flat, y_flat):
    # full log probability 
    lgl = logl(p, element, x, y, s, fmean_input, wlog_input, x_flat, y_flat)
    return np.sum(logprior(p, element) + lgl)


def runMCMC(element, wlog_input, fmean_input,
            x_flat, y_flat_sm, y_flat, y_flat_err,
            spec, posterior_save=False,
            plot_save=False, file_save=False, plotChain=False):

    '''runMCMC: sample probability distribution using package emcee,
    and get marginalized distribution of model parameters  '''

    x0 = X0[element]['cuts']
    p00 = P0[element]
    prior = Prior[element]

    ndim, nwalkers = 4, 6 * 2

    # template fit region
    lower_inds = (x_flat > x0[0]) & (x_flat < x0[1])
    lower_x_flat = x_flat[lower_inds]
    lower_y_flat_sm = y_flat_sm[lower_inds]
    lower = np.min(lower_x_flat[lower_y_flat_sm ==
                      np.max(lower_y_flat_sm)]) + 100
    upper_inds = (x_flat > x0[2]) & (x_flat < x0[3])
    upper_x_flat = x_flat[upper_inds]
    upper_y_flat_sm = y_flat_sm[upper_inds]
    upper = np.max(upper_x_flat[upper_y_flat_sm ==
                      np.max(upper_y_flat_sm)]) - 100

    elx = x_flat[(x_flat > lower) & (x_flat < upper)]
    ely = y_flat[(x_flat > lower) & (x_flat < upper)]
    els = y_flat_err[(x_flat > lower) & (x_flat < upper)]

    best_pos = []
    p0 = [p00 + 1e-6 * np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, logp,
                                      args=(element, elx, ely, els,
                                            fmean_input, wlog_input,
                                            x_flat, y_flat))
    # run MCMC for 30 steps starting from the tiny ball defined above
    sampler.run_mcmc(p0, 30)
    best_pos.append(sampler.flatchain[sampler.flatlnprobability.argmax()])
    pos = emcee.utils.sample_ball(best_pos[-1], best_pos[-1] / 100.,
                                  size=nwalkers)
    sampler.reset()
    sampler.run_mcmc(pos, 1000)
    best_pos.append(sampler.flatchain[sampler.flatlnprobability.argmax()])
    value_50 = [np.percentile(sampler.chain[:
        , :
        , 0], [50])[0],
                np.percentile(sampler.chain[:
        , :
        , 1], [50])[0],
                np.percentile(sampler.chain[:
        , :
        , 2], [50])[0],
                np.percentile(sampler.chain[:
        , :
        , 3], [50])[0]]

    # save marginalized distribution of model parameters
    if posterior_save:
        with open(outdir + "/" + posterior_save, 'w') as f:
            pkl.dump(sampler, f)

    # save template fit plot, corner plot (and chain plot if requested)
    if plot_save:

        # save template fit plot
        fig, ax = pl.subplots(figsize=(15, 15))
        
        fittemplate(value_50, element, fmean_input, wlog_input,
                    elx, ely, els, x_flat,
                    y_flat, ax = ax)
        
        pl.savefig(outdir + "/" + plot_save)
        pl.close(fig)

        # save corner plot
        value_50 = [np.percentile(sampler.chain[:
            , :
            , 0], [50])[0],
                    np.percentile(sampler.chain[:
            , :
            , 1], [50])[0] * 4,
                    np.percentile(sampler.chain[:
            , :
            , 2], [50])[0],
                    np.percentile(sampler.chain[:
            , :
            , 3], [50])[0]]

        sampler.flatchain[:
            , 1] = sampler.flatchain[:
            , 1] * 4

        rcParams['xtick.labelsize'] = 22.
        rcParams['ytick.labelsize'] = 22.
        fig_corner = corner.corner(sampler.flatchain, truths=value_50,
                      quantiles=[0.16, 0.5, 0.84],
                      labels=[r"$v$ [10$^3$ km s$^{-1}$]",
                              "$\sigma$ [$10^3$ km s$^{-1}$]", "$a$",
                              "$\Delta$$w$ [$\AA$]"])
        fig_corner.savefig(outdir + "/" + \
                           plot_save.replace('.pdf','Fit.pdf'))
        pl.close(fig_corner)

        # save chain plot
        if plotChain:
            y_label = [r"$v/1000$ [km s$^{-1}$]", "$\sigma/1000$ [km s$^{-1}$ ]",
                       "amplitude", "wave-range [\AA]"]
            figChain = pl.figure(figsize=(15, 4*ndim))
            for i in range(ndim):
                ax = figChain.add_subplot(ndim,1,i+1)
                ax.plot(range(1000), sampler.chain[:
                                                     , :
                                                     , i].T)
                if i == ndim - 1:
                    ax.set_xlabel("steps", fontsize=30)
                ax.set_ylabel(y_label[i], fontsize=30)
            figChain.savefig(outdir + "/" + \
                             plot_save.replace('.pdf','Chain.pdf'))
        pl.close('all')


    # save initial template fit region, mean acceptance fraction, initial
    # values for parameters, and 16th, 50th, 84th percentiles of marginalized
    # distribution of model parameters
    if file_save:
        f = open(outdir + "/" + file_save, 'w')
        f.write('region to find initial template fit region:' + str(x0) + '\n')
        f.write("Mean acceptance fraction: {0:.3f} \n"
                .format(np.mean(sampler.acceptance_fraction)))
        f.write('uniform prior for v/1000 in km/s, sigma/10 in angstrom, ' +
                'amplitude: ' + str(prior) + '\n')
        f.write('v/1000 in km/s, sigma/1000 in km/s, amplitude, ' +
                'wave-range in angstrom' + '\n')
        f.write('initial guess: ' + str(p00) + '\n')
        f.write('best value: ' + \
                ' '.join([str(bp) for bp in best_pos[-1]]) + '\n')
        f.write('16th, 50th, 84th percentiles \n')
        f.write(' '.join([str(pc) \
                          for pc in np.percentile(sampler.chain[:, :, 0],
                                                  [16, 50, 84])]) +
                ' for v/1000 in km/s\n')
        f.write(' '.join([str(pc) \
                          for pc in np.percentile(sampler.chain[:, :, 1],
                                                  [16, 50, 84])]) +
                         ' for sigma/1000 in km/s\n')
        f.write(' '.join([str(pc) \
                          for pc in np.percentile(sampler.chain[:, :, 2],
                                                  [16, 50, 84])]) +
                         ' for amplitude\n')
        f.write(' '.join([str(pc) \
                          for pc in np.percentile(sampler.chain[:, :, 3],
                                                  [16, 50, 84])]) +
                         ' for wave-range in angstrom\n')

        f.close()

    print('Mean acceptance fraction: {0:.3f}'
          .format(np.mean(sampler.acceptance_fraction)))
    print('{0:10} {1:15} {2:15} {3:5}  percentiles of marginalized distribution of model parameters'\
          .format(" ", "16th", "50th", "84th"))

    # 16th, 50th, 84th percentiles of the velocity/1000
    pargs = np.percentile(sampler.chain[:, :, 0], [16, 50, 84])
    print ('{0:15.3f} {1:15.3f} {2:15.3f}   for v/1000 in km/s'\
        .format(pargs[0], pargs[1], pargs[2]))
    
    # 16th, 50th, 84th percentiles of the sigma/10000 in km/s
    pargs = np.percentile(sampler.chain[:, :, 1], [16, 50, 84])
    print ('{0:15.3f} {1:15.3f} {2:15.3f}   for sigma/1000 in km/s'\
        .format(pargs[0], pargs[1], pargs[2]))

    # 16th, 50th, 84th percentiles of the amplitude
    pargs = np.percentile(sampler.chain[:, :, 2], [16, 50, 84])
    print ('{0:15.3f} {1:15.3f} {2:15.3f}   for amplitude'\
        .format(pargs[0], pargs[1], pargs[2]))

    # 16th, 50th, 84th percentiles of the wave-range in angstrom
    pargs = np.percentile(sampler.chain[:, :, 3], [16, 50, 84])
    print ('{0:15.3f} {1:15.3f} {2:15.3f}   for wavelenght range in A'\
        .format(pargs[0], pargs[1], pargs[2]))


def conv(spec, template, element):
    # fires off code to read in input files, fired off MCMC

    print ("Working with element: ", element)
    # check that the element has all necessary input info
    if element not in P0.keys() or element not in Prior.keys() or \
       element not in X0.keys():

        print("\nElement not enabled:")
        print("You must set up the initial parameters and ")
        print("prior parameters for %s in elementsDicts.py\n"%element)

        return -1
    
    
    # create output directory if it does not exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    try:
        wlog_input, fmean_input, x_flat, y_flat_sm, \
            y_flat, y_flat_err, phase = \
                                        readdata(spec, template)
    except IOError:
        print("\n\n You must pass a .csv and a .sav file, or 2 .sav ")
        print("files or a file and a phase as input")
        return -1
    if isinstance(wlog_input, int):
        return -1

    nameroot = spec.split('/')[-1].split('.')[0] + "-"  + element

    print("running convolution...")
    
    print('\n\nElement: {0}, Phase: {1:d}, input spectrum: {2}'\
          .format(element, phase, spec) )

    if y_flat[x_flat < X0[element]['min']].mean() and \
       y_flat[x_flat > X0[element]['max']].mean() :
        t1 = time.time()
        
        runMCMC(element, wlog_input, fmean_input, x_flat, y_flat_sm,
                y_flat, y_flat_err, spec, 
                posterior_save=nameroot + '.p',
                plot_save=nameroot + '.pdf',
                file_save=nameroot + '.dat',
                plotChain=False)
        t2 = time.time()
        print('minimization took {} seconds'.format(t2 - t1))

        
    else:
        print("\nError: wavelength range doesn't match")
        return -1
    return 0
if __name__ == '__main__':

    helpStrg = '''Arguments: observed spectrum and template spectrum: use as
$python Ic_conv_Icbl_MCMC.py 10qts_20100815_Lick_3-m_v1-z.flm-flat.sav  meanspecIc_0.sav Fe
or
$python Ic_conv_Icbl_MCMC.py 10qts_20100815_Lick_3-m_v1-z.flm-flat.sav  0 Fe'''

    element = 'Fe'
    if len(sys.argv) == 1:
        # use default arguments for testing
        print ("\n\n Hallo!\n\nThis is a test using SN PTF10qts at phase 0")
        obsSpec = '10qts_20100815_Lick_3-m_v1-z.flm-flat.sav'
        templSpec = 'IcTemplates/meanspecIc_0.sav'

    elif len(sys.argv) >= 3:
        # assume that you passed correctly the spectrum and template spectrum
        obsSpec, templSpec = sys.argv[1:3]
        
        if len(sys.argv) == 4:
            # element to fit
            element = sys.argv[3]                

    else:
        print(helpStrg)
        sys.exit()

    if conv(obsSpec, templSpec, element) == -1:
        print(helpStrg)
        sys.exit()
        
