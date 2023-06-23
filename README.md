[![DOI](https://zenodo.org/badge/22593/nyusngroup/SESNspectraLib.svg)](https://zenodo.org/badge/latestdoi/22593/nyusngroup/SESNspectraLib)

# SESNspectraLib

This repository contains code that was developed and used in [Liu et al. (2016)](http://adsabs.harvard.edu/abs/2016ApJ...827...90L) and [Modjaz et al. (2016)](http://adsabs.harvard.edu/abs/2016ApJ...832..108M) to process spectra of Stripped Envelope Supenova (SESNe), i.e., SNe of types IIb, Ib, Ic and Ic-bl. 
Common issues encountered in using SN spectral data include missing uncertainties and inhomogeneous measurements of the velocity for differen SN types. Our library intends to address these issues. 

At this time two pieces of code are provided: 
- **SNspecFFTsmooth.pro**: IDL routine to smooth SN spectra by separating SN signal from noise in Fourier space and to output the constructed spectral uncertainty array (AKA "error spectrum"). It should be used to calculate the uncertainty of spectra in absence of reduction-generated spectral uncertainty arrays, or to replace those in a comprehensive manner for a large and diverse dataset ([Liu et al. 2016](http://adsabs.harvard.edu/abs/2016ApJ...827...90L)). Note that galaxy emissions should be clipped beforehand. 
- **Ic_conv_Icbl_MCMC.py**: Python module to measure the absorption and width velocities of SN Ic-bl or SLSNe Ic from the Fe II 5169 spectral feature. This code provides and implements a method to measure the spectral velocity of SN Ic-bl consistently with what is usually done for other subtypes, for fair and meanungful comparisons, enabling the measurement of the Fe II 5169 blue-shift despite the fact that this feature is blended in Ic-bl specra (typically the Fe II triplet appears as a single blended feature in Ic-bl). The code measures the relative blue-shift and width of Fe II 5169 in SN Ic-bl spectrum *relative to the SN Ic mean template at the same phase*. Note that code outputs the width of the Gaussian broadening kernel as sigma, and thus has to be multiplied by 2.354 to be converted to a FWHM, as computed in Modjaz et al. (2016) for SNe Ic-bl and in Liu, Modjaz & Bianco (2017) for SLSNe Ic. This is done by matching the SN Ic-bl spectrum under consideration to a broadened and blueshifted SN Ic mean template at the same phase. The absorption velocities of the SNe Ic templates and their uncertainties are in [mean_FeII5169_vabs_comb.dat](https://github.com/nyusngroup/SESNspectraLib/blob/master/mean_FeII5169_vabs_comb.dat). In order to obtain *absolute absorption velocities* the outputted relative absorption velocity of the fit has to be *added* to the template Ic absorption velocity at the chosen phase, andn the uncertainties (one from the fit and the other from the standard deviation of all the measured velocities of SNe Ic that went into the Ic template - to represent the diversity of Ic vels at that phase) have to be added *in quadrature*. We used the code for measuring velocities for the feature Fe II 5169 ([Modjaz et al. 2016](http://adsabs.harvard.edu/abs/2016ApJ...832..108M)) but in principle, it can be used for other features as well. See [here](https://github.com/nyusngroup/SESNspectraLib/blob/master/M16AppendixA.pdf) (Appendix A of Modjaz et al. 2016) for a detailed discussion of this method.


NOTE: When not specified otherwise (e.g. Ic_conv_Icbl) our code can be used on any medium resolution spectra of SNe of any type, although it is designed for and tested on SESN and SLSNe Ic.

### IDL Libraries
The Meron IDL library can be found publicly available [here.](https://millenia.cars.aps.anl.gov/software/idl/index.html)

### Examples
1) The following chunk of IDL code shows how to run **SNspecFFTsmooth.pro** to produce the top panel of Figure 17 in [Liu et al. (2016)](http://adsabs.harvard.edu/abs/2015arXiv151008049L):
```
IDL> readcol, 'sn2004gq-20041212-z-bl.flm', w, f
IDL> SNspecFFTsmooth, w, f, 1000, f_ft, f_std, sep_vel
IDL> plot, w, f
IDL> oplot, w, f_ft, color=250
```
![alt tag](https://github.com/nyusngroup/SESNspectraLib/blob/master/example_IDL_plot.png)

2) The following chunk of code shows how to run **Ic_conv_Icbl_MCMC.py**.
In the following example, the input file is "10qts_20100815_Lick_3-m_v1-z.flm-flat.csv" (see header in the code for format) and the code produces the following outputs:
- **10qts_20100815_Lick_3-m_v1-z.flm-flat.sav-Fe.dat**, a text file that contains initial template fit region, **emcee** mean acceptance fraction, initial values for parameters, and 0.15th, 2.5th, 16th, 50th, 84th, 97.5th, 99.85th percentiles of marginalized distribution of model parameters.
- **10qts_20100815_Lick_3-m_v1-z.flm-flat.sav-Fe.p**, a pickle file that contains the marginalized distributions of model parameters.
- **10qts_20100815_Lick_3-m_v1-z.flm-flat.sav-Fe.pdf**, a PDF file that contains template fit plot
- **10qts_20100815_Lick_3-m_v1-z.flm-flat.sav-FeFit.pdf** corner plot (produced using the [corner.py] (https://github.com/dfm/corner.py package) if the package is installed) as a fit diagnostic.
- **10qts_20100815_Lick_3-m_v1-z.flm-flat.sav-FeChain.pdf**, optional, a chain plot showing the MCMC walkers as a fit diagnostic (http://dan.iel.fm/emcee/current/).

Note that the input spectra file need to contain the flattened spectrum which can be produced with **snidflat.pro** (using the default values) in **snid_pro.tgz**, downloadable at [Stephane Blondin's webpage](https://people.lam.fr/blondin.stephane/software/snid/index.html#Download), and the smoothed spectrum and uncertainty array, which can be produced with **SNspecFFTsmooth.pro**.

```
>>> import Ic_conv_Icbl_MCMC as Ic_conv_Icbl
>>> Ic_conv_Icbl.conv('10qts_20100815_Lick_3-m_v1-z.flm-flat.csv', 0, 'Fe')
````
or from the command line 
``` 
$python Ic_conv_Icbl_MCMC.py 10qts_20100815_Lick_3-m_v1-z.flm-flat.csv 0 Fe
Working with element:  Fe
reading inputs...
running convolution...


emcee Mman acceptance fraction: 0.576
           16th            50th            84th   percentiles of marginalized distribution of model parameters
         11.010          11.297          11.560   BLUE-SHIFT v/1000 in km/s WITH RESPECT to the Ic template at the chosen phase
          4.099           4.471           4.792   BROADENING sigma/1000 in km/s  WITH RESPECT to the Ic template at the chosen phase
          1.302           1.400           1.486   NORMALIZATION amplitude
          0.731           1.981           5.907   WAVELENGTH RANGE ADJUSTMENT from 470.19 in A
minimization took 9.92309498787 seconds

```
The main plot out will look like this:
![alt tag](https://raw.githubusercontent.com/nyusngroup/SESNspectraLib/master/10qts_20100815_Lick_3-m_v1-z.flm-flat-Fe.png)



###Required packages:
This code makes use of several default python packages, including [NumPy](http://www.numpy.org/), [SciPy](https://www.scipy.org/) (version >=0.10), [Matplotlib](http://matplotlib.org/) and [pickle](https://docs.python.org/2/library/pickle.html), and a few additional, less common, packages:

- [emcee](http://dan.iel.fm/emcee/current/)
- [corner](https://github.com/dfm/corner.py) (optional, but recommanded: it generates plots of the marginalized distributions as visual diagnostics of your fit)

###Acknowledgement:

If you use code in this repository, please <b>acknowledge</b> this work by including in your paper:

	This work made use of the data products generated by the NYU SN group, and 
	released under DOI:10.5281/zenodo.58767, 
	available at \url{https://github.com/nyusngroup/SESNspectraLib}.
	  
*and* <b>cite</b> the appropriate reference(s):


Modjaz et al. (2016):

    @ARTICLE{2016ApJ...832..108M,
   	author = {{Modjaz}, M. and {Liu}, Y.~Q. and {Bianco}, F.~B. and {Graur}, O.
	},
    	title = "{The Spectral SN-GRB Connection: Systematic Spectral Comparisons between Type Ic Supernovae and Broad-lined Type Ic Supernovae with and without Gamma-Ray Bursts}",
 	journal = {\apj},
	archivePrefix = "arXiv",
	eprint = {1509.07124},
	primaryClass = "astro-ph.HE",
	keywords = {gamma-ray burst: general, gamma-ray burst: individual: GRB-980425, GRB-030329, GRB-060218, GRB-100316D, GRB-120422A, GRB-	130427A, GRB-130702A, GRB-130215A, supernovae: general, supernovae: individual: SN-1994I, SN-2004aw, SN-2007gr, SN-1998bw, SN-2003dh, SN-2006aj, SN-2009bb, SN-2010bh, SN-2012ap, SN-2012bz, SN-2013cq, SN-2013dx, SN-2013ez },
     	year = 2016,
    	month = dec,
   	volume = 832,
      	eid = {108},
    	pages = {108},
      	doi = {10.3847/0004-637X/832/2/108},
  	adsurl = {http://adsabs.harvard.edu/abs/2016ApJ...832..108M},
 	adsnote = {Provided by the SAO/NASA Astrophysics Data System}
	}
    

Liu et al. (2016):

    @ARTICLE{2016ApJ...827...90L,
   	author = {{Liu}, Y.-Q. and {Modjaz}, M. and {Bianco}, F.~B. and {Graur}, O.
	},
    	title = "{Analyzing the Largest Spectroscopic Data Set of Stripped Supernovae to Improve Their Identifications and Constrain Their Progenitors}",
  	journal = {\apj},
	archivePrefix = "arXiv",
   	eprint = {1510.08049},
 	primaryClass = "astro-ph.HE",
 	keywords = {methods: data analysis, supernovae: general, supernovae: individual: SNe 1993J, 2005bf, 2005E, 2011dh},
     	year = 2016,
    	month = aug,
   	volume = 827,
      	eid = {90},
    	pages = {90},
      	doi = {10.3847/0004-637X/827/2/90},
   	adsurl = {http://adsabs.harvard.edu/abs/2016ApJ...827...90L},
  	adsnote = {Provided by the SAO/NASA Astrophysics Data System}
	}
