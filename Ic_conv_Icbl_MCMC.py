# -*- coding: utf-8 -*-
# readin: flattened Ic-bl spectrum and the corresponding uncertainty array, SNe Ic template
# output: marginalized distribution of model parameters, including but not limited to absorption velocity
# note: prior and initial values of model parameters, and initial template fitting region can be changed as needed

import sys
import os
import numpy as np
from scipy.io.idl import readsav
import pylab as pl
from scipy.ndimage import filters
from scipy.signal import  gaussian
from scipy.interpolate import interp1d
from scipy.optimize import minimize 
import time
from matplotlib.backends.backend_pdf import PdfPages
import triangle 
import emcee
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import pickle as pkl

t1 = time.time()



# read in flattened Ic-bl spectrum and the corresponding uncertainty array, SNe Ic template
def readdata(spec, template):

    # read in Ic template            
    s=readsav(template) # read in Ic template
    wlog_input=s.wlog[np.where((s.wlog > 4400) & (s.wlog < 9000))]
    fmean_input=s.fmean[np.where((s.wlog > 4400) & (s.wlog < 9000))]
    
    # read in Icbl spectrum and uncertainty array
    s2=readsav(spec) # read in Icbl spectrum
    x_flat=s2.wavelog_input[0:1024]       
    y_flat_sm=s2.flatflux_input_sm
    y_flat=s2.flatflux_input
    y_flat_err=s2.flatflux_err_input

    return wlog_input,fmean_input, x_flat,y_flat_sm,y_flat,y_flat_err   
    
# fit SN Ic-bl spectrum to a convolved and blueshifted SN Ic template
def fittemplate(p,fmean_input,wlog_input,lx,ly,ly_err, plot=False):
    stre,scale1,v,sig=p[3], p[2]*100,-p[1]*1000,p[0]*10
    w_lower=lx[0]+scale1
    w_upper=lx[-1]-scale1

    inds=(lx>w_lower) * (lx<w_upper)
    ly_new=ly[inds]
    ly_err_new=ly_err[inds]
    lx_new=lx[inds]

    beta=v/299792.458
    doppler=np.sqrt((1+beta)/(1-beta))
    b = gaussian(300, sig)
    thisy=filters.convolve1d(stre*fmean_input, b/b.sum())  # first do Gaussian convolution
    f2 = interp1d(wlog_input*doppler, thisy,bounds_error=False, fill_value=0)(lx_new) 
    chisq=np.sum((ly_new-f2)**2/ly_err_new**2)/(len(ly_new)-len(p))
    
    if plot:
        pl.figure(figsize=(15,15))
        #pl.plot(x_flat,y_flat,'k',alpha=0.5,label=r"$t=$%.2f days"%(float(0.0)))
        pl.plot(lx_new,ly_new,'k')
        pl.plot(lx_new,f2,'r',linewidth=3)
        pl.plot(lx[0],ly[0],'o',color='blue',)
        pl.plot(lx[-1],ly[-1],'o',color='blue',)
        pl.plot(wlog_input*doppler, stre*fmean_input,'r',linewidth=2,alpha=0.5,label=r"$\sigma$=%.2f \AA~  v=%.2f km/s~  w=%.2f \AA~  a=%.2f~ $\chi^2_r$=%.2f"%(sig,-v,scale1,stre,chisq))
        pl.xlabel("wavelength (\AA)",fontsize=25)
        pl.ylabel("flux (normalized)",fontsize=25)
        pl.legend(fontsize=25)
            
    return chisq
    
# log prior        
def logprior (p):
    s=p[0]
    v=p[1]
    scale1=p[2]
    stre=p[3]
    #scale2=p[3]
    if s>priorFe[0] and s<priorFe[1] and v>priorFe[2] and v<priorFe[3] and  stre>priorFe[6] and stre<priorFe[7]:# and scale1>priorFe[4] and scale1<priorFe[5]:
        #return 0.0+np.log(1/(0.33 * np.sqrt(2 * np.pi)) * np.exp( - (scale1 - 1)**2 / (2 * 0.33**2) ))
        return np.log(np.exp( - (scale1 - 0.1)**2 / (2 * 0.5**2) ))
    return -np.inf

# log likelihood                    
def logl(p,x,y,s,fmean_input,wlog_input):
    return  -np.log(s)-0.5*(fittemplate(p,fmean_input,wlog_input,x,y,s))
    
# full log probability     
def logp(p,x,y,s,fmean_input,wlog_input):
    lgl= logl(p,x,y,s,fmean_input,wlog_input)
    #print lgl
    return np.sum(logprior(p)+lgl)

# sample probability distribution using package emcee, and get marginalized distribution of model parameters  
def runMCMC(wlog_input,fmean_input, x_flat,y_flat_sm,y_flat,y_flat_err,x0,p00,prior,spec,posterior_save='no',plot_save='no',file_save='no'):
    ndim, nwalkers = 4, 6*2
    
    Fe_lower_inds=(x_flat>x0[0]) * (x_flat<x0[1])
    Fe_lower_x_flat=x_flat[Fe_lower_inds]
    Fe_lower_y_flat_sm=y_flat_sm[Fe_lower_inds]
    #Fe_lower_y_flat=y_flat[Fe_lower_inds]                        
    Fe_lower=np.min(Fe_lower_x_flat[np.where(Fe_lower_y_flat_sm == np.max(Fe_lower_y_flat_sm))])+100
    Fe_upper_inds=(x_flat>x0[2]) * (x_flat<x0[3])
    Fe_upper_x_flat=x_flat[Fe_upper_inds]
    Fe_upper_y_flat_sm=y_flat_sm[Fe_upper_inds]
    #Fe_upper_y_flat=y_flat[Fe_upper_inds]
    Fe_upper=np.max(Fe_upper_x_flat[np.where(Fe_upper_y_flat_sm == np.max(Fe_upper_y_flat_sm))])-100
                        
    Fex=x_flat[np.where((x_flat>Fe_lower) & (x_flat<Fe_upper))]
    Fey=y_flat[np.where((x_flat>Fe_lower) & (x_flat<Fe_upper))]
    Fey_err=y_flat_err[np.where((x_flat>Fe_lower) & (x_flat<Fe_upper))]
                                            
    Fes=Fey_err #np.ones(len(Fex))/np.sqrt(len(Fex))

    best_pos=[]
    p0 = [p00 + 1e-6*np.random.randn(ndim) for i in range(nwalkers)]
    samplerFe = emcee.EnsembleSampler(nwalkers, ndim, logp, 
                                            args=(Fex,Fey,Fes,fmean_input,wlog_input))
    samplerFe.run_mcmc(p0, 30) # run MCMC for 30 steps starting from the tiny ball defined above
    best_pos.append(samplerFe.flatchain[samplerFe.flatlnprobability.argmax()])
    pos = emcee.utils.sample_ball(best_pos[-1], best_pos[-1]/100., size=nwalkers)
    samplerFe.reset()
    samplerFe.run_mcmc(pos, 1000)
    best_pos.append(samplerFe.flatchain[samplerFe.flatlnprobability.argmax()])
    value_50=[np.percentile(samplerFe.chain[:,:,0],[50])[0],np.percentile(samplerFe.chain[:,:,1],[50])[0], # median values
            np.percentile(samplerFe.chain[:,:,2],[50])[0],np.percentile(samplerFe.chain[:,:,3],[50])[0]]

    # save marginalized distribution of model parameters
    if posterior_save !='no':
        with open(posterior_save,'w') as f:
            pkl.dump(samplerFe,f)
    
    # save corner plots
    if plot_save != 'no':  
        pp = PdfPages(plot_save)
                
        #fittemplate(best_pos[-1],fmean_input,wlog_input,Fex,Fey,Fes,plot=True)
        fittemplate(value_50,fmean_input,wlog_input,Fex,Fey,Fes,plot=True)
        pp.savefig()
        pl.close() # very important to close figures otherwise may comsume too much memory

        #samplerFe = pkl.load( open(posterior_save , "rb" ) ) #this line shows how to use the stored pickle file to plot posterior distribution        
        triangle.corner(samplerFe.flatchain,truths=value_50,labels=[r"$\sigma$ [10 \AA]", "$v$ [10$^3$ km s$^{-1}$]", "$w$ [10$^2$ \AA]", "$a$"])    
        pp.savefig()
        pl.close() # very important to close figures otherwise may comsume too much memory

        pp.close()

    # save initial template fit region, mean acceptance fraction, initial values for parameters and others
    if file_save != 'no':   
        f=open(file_save,'w')        
        f.write('region: '+str(x0)+' to find initial template fit region: '+str(np.array([round(Fe_lower),round(Fe_upper)]))+'\n')
        f.write("Mean acceptance fraction: {0:.3f} \n"
                        .format(np.mean(samplerFe.acceptance_fraction)))    
        f.write('sigma/10 v/1000 scale/100 stretch, prior: '+str(prior)+'\n')
        f.write('initial guess: '+str(p00)+'\n')
        f.write('best value: '+str(best_pos[-1])+'\n')
        f.write('16th, 50th, 84th percentiles \n')
        f.write(str(np.percentile(samplerFe.chain[:,:,0],[16,50,84]))+'\n')
        f.write(str(np.percentile(samplerFe.chain[:,:,1],[16,50,84]))+'\n')
        f.write(str(np.percentile(samplerFe.chain[:,:,2],[16,50,84]))+'\n')
        f.write(str(np.percentile(samplerFe.chain[:,:,3],[16,50,84]))+'\n')
        f.close()                   

    print("Mean acceptance fraction: {0:.3f}"
                    .format(np.mean(samplerFe.acceptance_fraction)))    
    print best_pos[-1] # sigma/10, velocity/1000, scale/100, stretch
    print np.percentile(samplerFe.chain[:,:,0],[16,50,84]) # 16th, 50th, 84th percentiles of the sigma/10 => 1 sigma error bar
    print np.percentile(samplerFe.chain[:,:,1],[16,50,84]) # 16th, 50th, 84th percentiles of the velocity/1000
    print np.percentile(samplerFe.chain[:,:,2],[16,50,84]) # 16th, 50th, 84th percentiles of the scale/100
    print np.percentile(samplerFe.chain[:,:,3],[16,50,84]) # 16th, 50th, 84th percentiles of the stretch        


def conv(spec, template):
                                                                                
    wlog_input,fmean_input, x_flat,y_flat_sm,y_flat,y_flat_err=readdata(spec, template)
    if np.mean(y_flat[np.where(x_flat < 4500)]) !=0 and np.mean(y_flat[np.where(x_flat > 5100)]) !=0: 
        runMCMC(wlog_input,fmean_input, x_flat,y_flat_sm,y_flat,y_flat_err,x0Fe,p0Fe,priorFe,spec,posterior_save=spec+'-Fe.p',plot_save=spec+'-Fe.pdf',file_save=spec+'-Fe.dat')                      
    else:
        print("wavelength range doesn't match")                

# initial parameter values
p0Fe=np.array([1,11,0.1,1.5]) #sigma/10, v/1000, wave-range/100, y-amplitude
# region to find initial template fit region
x0Fe=np.array([4200,4800,5000,5600])                    
# prior
priorFe=np.array([0,5,0,40,0,4,0,3]) # to remove 0 and 4

conv('10qts_20100815_Lick_3-m_v1-z.flm-flat.sav', 'meanspecIc_0.sav')

#spec = '10qts_20100815_Lick_3-m_v1-z.flm-flat.sav'
#template = 'meanspecIc_0.sav'
#wlog_input,fmean_input, x_flat,y_flat_sm,y_flat,y_flat_err=readdata(spec, template)
#if np.mean(y_flat[np.where(x_flat < 4500)]) !=0 and np.mean(y_flat[np.where(x_flat > 5100)]) !=0: 
#    runMCMC(wlog_input,fmean_input, x_flat,y_flat_sm,y_flat,y_flat_err,x0Fe,p0Fe,priorFe,spec,posterior_save=spec+'-Fe.p',plot_save=spec+'-Fe.pdf',file_save=spec+'-Fe.dat')                      
#else:
#    print("wavelength range doesn't match")     
#        
t2 = time.time()
print 'minimization took {} seconds'.format(t2 - t1)

