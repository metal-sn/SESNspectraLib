# -*- coding: utf-8 -*-
# readin: flattened Ic-bl spectrum and the corresponding uncertainty array, SNe Ic template
# output: marginalized distribution of model parameters, including but not limited to absorption velocity, see Modjaz et al. (2016) for details
# note: prior and initial values of model parameters, and initial template fitting region can be changed as needed

import numpy as np
from scipy.io.idl import readsav
import pylab as pl
from scipy.ndimage import filters
from scipy.signal import  gaussian
from scipy.interpolate import interp1d
from matplotlib.backends.backend_pdf import PdfPages
import corner 
import emcee
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl
import pickle as pkl
import time

# initial parameter values
#v/1000, sigma/10, y-amplitude, wave-range/10
P0Fe=np.array([11.0,1.0,1.0,1.0]) 
# region to find initial template fit region
X0Fe = np.array([4200,4800,5000,5600])                    
# prior for v/1000, sigma/10, y-amplitude
PriorFe = np.array([0,30,0,5,0,3])
    
def readdata(spec, template):
    ''' read in flattened Ic-bl spectrum and the corresponding 
    uncertainty array,  SNe Ic template '''

    # read in Ic template            
    s = readsav(template)
    wlog_input = s.wlog[np.where((s.wlog > 4400) & (s.wlog < 9000))]
    fmean_input = s.fmean[np.where((s.wlog > 4400) & (s.wlog < 9000))]
    
    # read in Icbl spectrum and uncertainty array
    s2 = readsav(spec)
    x_flat = s2.wavelog_input[0:1024]       
    y_flat_sm = s2.flatflux_input_sm
    y_flat = s2.flatflux_input
    y_flat_err = s2.flatflux_err_input

    return wlog_input,fmean_input, x_flat,y_flat_sm,y_flat,y_flat_err   
    
def fittemplate(p,fmean_input,wlog_input,lx,ly,ly_err, x_flat, y_flat, plot = False):
    ''' fit SN Ic-bl spectrum to a convolved and blueshifted SN Ic template'''
    v, sig, amplitude, w_range = -p[0]*1000, p[1]*10, p[2], p[3]*10
    w_lower = lx[0]+w_range
    w_upper = lx[-1]-w_range

    inds = (lx>w_lower) * (lx<w_upper)
    ly_new = ly[inds]
    ly_err_new = ly_err[inds]
    lx_new = lx[inds]

    beta = v/299792.458
    doppler = np.sqrt((1+beta)/(1-beta))
    b  =  gaussian(300, sig)

    # first do Gaussian convolution    
    thisy = filters.convolve1d(amplitude*fmean_input, b/b.sum())  
    f2  =  interp1d(wlog_input*doppler, thisy,bounds_error = False, fill_value = 0)(lx_new) 
    chisq = np.sum((ly_new-f2)**2/ly_err_new**2)/(len(ly_new)-len(p))
    
    if plot:
        fig,ax = pl.subplots(figsize = (15,15))
        minorLocatory   = MultipleLocator(0.02)
        minorLocatorx   = MultipleLocator(100)
        ax.xaxis.set_minor_locator(minorLocatorx)
        ax.yaxis.set_minor_locator(minorLocatory)

        pl.plot(x_flat,y_flat,'k',alpha = 0.5)
        pl.plot(lx_new,ly_new,'k')
        pl.plot(lx_new,f2,'r',linewidth = 3)
        pl.plot(lx[0],ly[0],'o',color = 'blue',)
        pl.plot(lx[-1],ly[-1],'o',color = 'blue',)
        pl.plot(wlog_input*doppler, amplitude*fmean_input, 'r', linewidth = 2,
                alpha = 0.5)
        pl.text(2200,0.5,r"$v$=%.0f km~s$^{-1}$, $\sigma$=%.0f km s$^{-1}$, $a$=%.1f, $\Delta$$w$=%.0f \AA"%(-v,sig*400,amplitude,w_range), fontsize=35)
        pl.text(5500,0.3,r"$\chi^2_r$=%.1f"%(chisq), fontsize=35)
        pl.xlabel("Rest Wavelength (\AA)", fontsize=35)
        pl.ylabel("Relative Flux", fontsize=35)
        pl.legend(fontsize = 35)
            
    return chisq
    
def logprior (p):
    # log prior        
    v=p[0]
    s=p[1]
    amplitude=p[2]
    w_range=p[3]
    if s > PriorFe[2] and s < PriorFe[3] and v > PriorFe[0] and \
       v < PriorFe[1] and  amplitude > PriorFe[4] and amplitude < PriorFe[5]:
        return np.log(np.exp( - (w_range)**2 / (2 * 3.3**2) ))
    return -np.inf

def logl(p, x, y, s, fmean_input, wlog_input, x_flat, y_flat):
    # log likelihood                    
    return  -np.log(s) - 0.5 * (fittemplate(p, fmean_input, wlog_input,
                                            x, y, s, x_flat, y_flat))
    
def logp(p, x, y, s, fmean_input, wlog_input, x_flat, y_flat):
    '''full log probability '''    
    lgl= logl(p, x, y, s, fmean_input, wlog_input, x_flat, y_flat)
    return np.sum(logprior(p) + lgl)

def runMCMC(wlog_input, fmean_input, x_flat, y_flat_sm, y_flat,
            y_flat_err, spec, 
            x0 = X0Fe, p00 = P0Fe, prior=PriorFe,
            posterior_save = False,
            plot_save = False, file_save = False):
    ''' sample probability distribution using package emcee, 
    and get marginalized distribution of model parameters  '''

    ndim, nwalkers = 4, 6 * 2
    
    Fe_lower_inds = (x_flat > x0[0]) * (x_flat < x0[1])
    Fe_lower_x_flat = x_flat[Fe_lower_inds]
    Fe_lower_y_flat_sm = y_flat_sm[Fe_lower_inds]
    #Fe_lower_y_flat=y_flat[Fe_lower_inds]                        
    Fe_lower = np.min(Fe_lower_x_flat[np.where(Fe_lower_y_flat_sm == np.max(Fe_lower_y_flat_sm))]) + 100
    Fe_upper_inds = (x_flat > x0[2]) * (x_flat < x0[3])
    Fe_upper_x_flat = x_flat[Fe_upper_inds]
    Fe_upper_y_flat_sm = y_flat_sm[Fe_upper_inds]
    #Fe_upper_y_flat=y_flat[Fe_upper_inds]
    Fe_upper = np.max(Fe_upper_x_flat[np.where(Fe_upper_y_flat_sm == np.max(Fe_upper_y_flat_sm))]) - 100
                        
    Fex = x_flat[np.where((x_flat > Fe_lower) & (x_flat < Fe_upper))]
    Fey = y_flat[np.where((x_flat > Fe_lower) & (x_flat < Fe_upper))]
    Fey_err = y_flat_err[np.where((x_flat > Fe_lower) & (x_flat < Fe_upper))]
                                            
    Fes = Fey_err #np.ones(len(Fex))/np.sqrt(len(Fex))

    best_pos = []
    p0 = [p00 + 1e-6*np.random.randn(ndim) for i in range(nwalkers)]
    samplerFe = emcee.EnsembleSampler(nwalkers, ndim, logp, 
                                            args=(Fex, Fey, Fes,
                                                  fmean_input, wlog_input, x_flat, y_flat))
    # run MCMC for 30 steps starting from the tiny ball defined above
    samplerFe.run_mcmc(p0, 30) 
    best_pos.append(samplerFe.flatchain[samplerFe.flatlnprobability.argmax()])
    pos = emcee.utils.sample_ball(best_pos[-1], best_pos[-1]/100.,
                                  size = nwalkers)
    samplerFe.reset()
    samplerFe.run_mcmc(pos, 1000)
    best_pos.append(samplerFe.flatchain[samplerFe.flatlnprobability.argmax()])
    value_50 = [np.percentile(samplerFe.chain[:,:,0],[50])[0],
                np.percentile(samplerFe.chain[:,:,1],[50])[0], # median values
                np.percentile(samplerFe.chain[:,:,2],[50])[0],
                np.percentile(samplerFe.chain[:,:,3],[50])[0]]

    # save marginalized distribution of model parameters
    if posterior_save:
        with open(posterior_save,'w') as f:
            pkl.dump(samplerFe,f)
    
    # save corner plots
    if plot_save:  
        pp = PdfPages(plot_save)
                
        #fittemplate(best_pos[-1],fmean_input,wlog_input,Fex,Fey,Fes,plot=True)
        fittemplate(value_50, fmean_input, wlog_input, Fex, Fey, Fes, x_flat, y_flat, plot=True)
        pp.savefig()
        pl.close() # very important to close figures otherwise may comsume too much memory

        value_50 = [np.percentile(samplerFe.chain[:,:,0],[50])[0],
                np.percentile(samplerFe.chain[:,:,1],[50])[0]*4, # median values
                np.percentile(samplerFe.chain[:,:,2],[50])[0],
                np.percentile(samplerFe.chain[:,:,3],[50])[0]]

        samplerFe.flatchain[:,1]=samplerFe.flatchain[:,1]*4

        mpl.rcParams['xtick.labelsize'] = 22.
        mpl.rcParams['ytick.labelsize'] = 22.

        #samplerFe = pkl.load( open(posterior_save , "rb" ) ) #this line shows how to use the stored pickle file to plot posterior distribution        
        corner.corner(samplerFe.flatchain, truths=value_50, quantiles=[0.16, 0.5, 0.84],
                      labels = [r"$v$ [10$^3$ km s$^{-1}$]", "$\sigma$ [$10^3$ km s$^{-1}$]", "$a$","$\Delta$$w$ [10 \AA]"])    
        pp.savefig()
        
        # very important to close figures otherwise may comsume too much memory
        pl.close() 
        
        #fig,ax = pl.subplots(figsize = (15,15))
        y_label=[r"$v/1000$ [km s$^{-1}$]", "$\sigma/1000$ [km s$^{-1}$]", "amplitude", "wave-range/10 [\AA]"]
        for i in range(ndim):
            pl.figure(figsize=(15,4))
            pl.plot(range(1000),samplerFe.chain[:,:,i].T) # nwalkers, iterations, ndim
            pl.xlabel("steps", fontsize=30)
            pl.ylabel(y_label[i], fontsize=30)
            pp.savefig()
        pl.close() 

        pp.close()

    # save initial template fit region, mean acceptance fraction, initial values for parameters and others
    if file_save:   
        f=open(file_save, 'w')        
        f.write('region: '+str(x0)+' to find initial template fit region: '+str(np.array([round(Fe_lower),round(Fe_upper)]))+'\n')
        f.write("Mean acceptance fraction: {0:.3f} \n"
                        .format(np.mean(samplerFe.acceptance_fraction)))    
        f.write('v/1000 sigma/1000 amplitude wave-range/10, prior: '+str(prior)+'\n')
        f.write('initial guess: '+str(p00)+'\n')
        f.write('best value: '+str(best_pos[-1])+'\n')
        f.write('16th, 50th, 84th percentiles \n')
        f.write(str(np.percentile(samplerFe.chain[:,:,0], [16,50,84]))+'\n')
        f.write(str(np.percentile(samplerFe.chain[:,:,1], [16,50,84]))+'\n')
        f.write(str(np.percentile(samplerFe.chain[:,:,2], [16,50,84]))+'\n')
        f.write(str(np.percentile(samplerFe.chain[:,:,3], [16,50,84]))+'\n')
        f.close()                   

    print("Mean acceptance fraction: {0:.3f}"
                    .format(np.mean(samplerFe.acceptance_fraction)))    
    print best_pos[-1] # velocity/1000 (km/s), sigma/1000 (km/s), amplitude, wave-range/10
    # 16th, 50th, 84th percentiles of the velocity/1000
    print np.percentile(samplerFe.chain[:,:,0],[16,50,84])*1000 
    # 16th, 50th, 84th percentiles of the sigma/10 => 1 sigma error bar
    print np.percentile(samplerFe.chain[:,:,1],[16,50,84])*1000 
    # 16th, 50th, 84th percentiles of the amplitude
    print np.percentile(samplerFe.chain[:,:,2],[16,50,84]) 
    # 16th, 50th, 84th percentiles of the wave-range/10
    print np.percentile(samplerFe.chain[:,:,3],[16,50,84]) *10


def conv(spec, template):
                                                                                
    wlog_input,fmean_input, x_flat,y_flat_sm,y_flat,y_flat_err=readdata(spec, template)
    if np.mean(y_flat[np.where(x_flat < 4500)]) != 0 and np.mean(y_flat[np.where(x_flat > 5100)]) !=0: 
        runMCMC(wlog_input, fmean_input, x_flat, y_flat_sm,
                y_flat, y_flat_err, spec,
                #x0Fe = X0Fe, p0Fe = P0Fe, priorFe = PriorFe,
                posterior_save = spec + '-Fe.p', plot_save=spec+'-Fe.pdf',
                file_save = spec+'-Fe.dat')                      
    else:
        print("wavelength range doesn't match")                


if __name__ == '__main__':
        
    t1 = time.time()

    conv('10qts_20100815_Lick_3-m_v1-z.flm-flat.sav', 'meanspecIc_0.sav')
    
    t2 = time.time()
    print 'minimization took {} seconds'.format(t2 - t1)

