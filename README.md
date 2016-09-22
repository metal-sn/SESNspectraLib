[![DOI](https://zenodo.org/badge/22593/nyusngroup/SESNspectraLib.svg)](https://zenodo.org/badge/latestdoi/22593/nyusngroup/SESNspectraLib)

# SESNspectraLib

This repository contains code that was developed and used in [Liu et al. (2016)](http://adsabs.harvard.edu/abs/2015arXiv151008049L) and [Modjaz et al. (2016)](http://adsabs.harvard.edu/abs/2015arXiv150907124M) to process spectra of Stripped Envelope Supenova (SESNe), i.e., SNe of types IIb, Ib, Ic and Ic-bl. 
Common issues encountered in using SN spectral data include missing uncertainties and inhomogeneous measurements of the velocity for differen SN types. Our library intends to address these issues. 

At this time two pieces of code are provided: 
- **SNspecFFTsmooth.pro**: IDL routine to smooth SN spectra by separating SN signal from noise in Fourier space and to output the constructed spectral uncertainty array (AKA "error spectrum"). It should be used to calculate the uncertainty of spectra in absence of reduction-generated spectral uncertainty arrays, or to replace those in a comprehensive manner for a large and diverse dataset ([Liu et al. (2016)](http://adsabs.harvard.edu/abs/2015arXiv151008049L)).
- **Ic_conv_Icbl_MCMC.py**: Python module to measure the absorption velocity of SN Ic-bl from the Fe II 5169 spectral feature. This code provides and implements a method to measure the spectral velocity of SN Ic-bl consistently with what is usually done for other subtypes, for fair and meanungful comparisons, enabling the measurement of the Fe II 5169 blue-shift despite the fact that this feature is blended in Ic-bl specra (typically the Fe II triplet appears as a single blended feature in Ic-bl). The code measures the blue-shift and width of Fe II 5169 by matching the SN Ic-bl spectrum under consideration to a broadened and blueshifted SN Ic mean template at the same phase. We use it for Fe II 5169 ([Modjaz et al. (2016)](http://adsabs.harvard.edu/abs/2015arXiv150907124M)) but it can be used in principle for other features as well. See [here](https://github.com/nyusngroup/SESNspectraLib/blob/master/M16AppendixA.pdf) (Appendix A of Modjaz et al. 2016) for a detailed discussion of this method.


NOTE: When not specified otherwise (e.g. Ic_conv_Icbl) our code can be used on any medium resolution spectra of SNe of any type, although it is designed for and tested on SESN.

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
         11.010          11.297          11.560   BLUE-SHIFT v/1000 in km/s
          4.099           4.471           4.792   BROADENING sigma/1000 in km/s
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

    @ARTICLE{2015arXiv150907124M,
      author = {{Modjaz}, M. and {Liu}, Y.~Q. and {Bianco}, F.~B. and {Graur}, O.
	    },
        title = "{The Spectral SN-GRB Connection: Systematic Spectral Comparisons between Type Ic Supernovae, and broad-lined Type Ic Supernovae with and without Gamma-Ray Bursts}",
      journal = {ArXiv e-prints},
    archivePrefix = "arXiv",
      eprint = {1509.07124},
    primaryClass = "astro-ph.HE",
    keywords = {Astrophysics - High Energy Astrophysical Phenomena, Astrophysics - Solar and Stellar Astrophysics},
         year = 2016,
        month = sep,
      adsurl = {http://adsabs.harvard.edu/abs/2015arXiv150907124M},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System} 
    }

Liu et al. (2016):

	@ARTICLE{2015arXiv151008049L,
   		author = {{Liu}, Y.-Q. and {Modjaz}, M. and {Bianco}, F.~B. and {Graur}, O.
		},
    		title = "{Analyzing the Largest Spectroscopic Dataset of Stripped Supernovae to Improve Their Identifications and Constrain Their Progenitors}",
  		journal = {ArXiv e-prints},
		archivePrefix = "arXiv",
		   eprint = {1510.08049},
 		primaryClass = "astro-ph.HE",
 		keywords = {Astrophysics - High Energy Astrophysical Phenomena, Astrophysics - Solar and Stellar Astrophysics},
     		year = 2016,
    		month = oct,
   		adsurl = {http://adsabs.harvard.edu/abs/2015arXiv151008049L},
  		adsnote = {Provided by the SAO/NASA Astrophysics Data System}
		}
