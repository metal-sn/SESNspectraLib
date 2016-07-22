# SESNspectraLib

This repository contains code that was used in [Liu et al. (2016)](http://adsabs.harvard.edu/abs/2015arXiv151008049L) and [Modjaz et al. (2016)](http://adsabs.harvard.edu/abs/2015arXiv150907124M) to process spectra of Stripped Envelope Supenova (SESNe),i.e., SNe of types IIb, Ib, Ic and Ic-bl. 
Common issues encountered in using SN spectral data include missing uncertainties and inhomogeneous measurements of the velocity for differen SN types. Our library intends to address these issues. 

At this time two pieces of code are provided: 
- **SNspecFFTsmooth.pro**: IDL routine to smooth SN spectra by separating SN signal from noise in Fourier space. It should be used to evaluate the uncertainty of spectra in absence of reduction-generated spectral uncertainty arrays, or to replace those in a comprehensive manner for a large and diverse dataset ([Liu et al. (2016)](http://adsabs.harvard.edu/abs/2015arXiv151008049L)).
- **Ic_conv_Icbl_MCMC.py**: Python module to measure the absorption velocity of SN Ic-bl from the Fe II 5169 spectral feature. This code provides and implements a method to measure the spectral velocity of SN Ic-bl consistently with what is usually done for oher subtypes, for fair and meanungful comparisons, enabling the measurement of the Fe II 5169 blue-shift despite the fact that this feature is blended in Ic-bl specra (typically the Fe II triplet appears as a single blended feature in Ic-bl). The code measures the blue-shift and width of Fe II 5169 by matching the SN Ic-bl spectrum under consideration to a broadened and blueshifted SN Ic mean template at similar phases. We use it for Fe II 5169 ([Modjaz et al. (2016)](http://adsabs.harvard.edu/abs/2015arXiv150907124M)) but it can be used in principle for other features as well.


NOTE: When not specified otherwise (e.g. Ic_conv_Icbl) our code can be used on any medium resolution spectra of SNe of any type, although it is designed for and tested on SESN.

### Examples
1) The following chunk of IDL code is an example of **SNspecFFTsmooth.pro** in IDL, which will produce the top panel of Figure 17 in [Liu et al. (2016)](http://adsabs.harvard.edu/abs/2015arXiv151008049L):
```
IDL> readcol, 'sn2004gq-20041212-z-bl.flm', w, f
IDL> SNspecFFTsmooth, w, f, 1000, w_ft, f_ft, f_std, sep_vel
IDL> plot, w, f
IDL> oplot, w_ft, f_ft, color=250
```
![alt tag](https://github.com/nyusngroup/SESNspectraLib/blob/master/example_IDL_plot.png)

2) The following chunk of Python is an example of **Ic_conv_Icbl_MCMC.py** in Python.
Note that the input spectra of **Ic_conv_Icbl_MCMC.py** are flattened by **snidflat.pro** (using the default values) in **snid_pro.tgz** which can be downloaded via [Stephane Blondin's webpage](https://people.lam.fr/blondin.stephane/software/snid/index.html#Download).

In the following example, the input is "10qts_20100815_Lick_3-m_v1-z.flm-flat.csv" (see header in the code for format) and the code produces the following outputs:
- **10qts_20100815_Lick_3-m_v1-z.flm-flat.sav-Fe.dat**, a text file that contains initial template fit region, mean acceptance fraction, initial values for parameters, and 16th, 50th, 84th percentiles of marginalized distribution of model parameters.
- **10qts_20100815_Lick_3-m_v1-z.flm-flat.sav-Fe.p**, a pickle file that contains marginalized distribution of model parameters.
- **10qts_20100815_Lick_3-m_v1-z.flm-flat.sav-Fe.pdf**, a PDF file that contains template fit plot
- **10qts_20100815_Lick_3-m_v1-z.flm-flat.sav-FeFit.pdf** corner plot (produced using the [corner.py] (https://github.com/dfm/corner.py package) if the package is installed) as a fit diagnosis
- **10qts_20100815_Lick_3-m_v1-z.flm-flat.sav-FeChain.pdf** as an optional plot, and chain plot for the MCMC fit, as a fit diagnosis.


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


Mean acceptance fraction: 0.560
           16th            50th            84th   percentiles of marginalized distribution of model parameters
         10.989          11.265          11.579   for v/1000 in km/s
          4.087           4.467           4.810   for sigma/1000 in km/s
          1.290           1.391           1.496   for amplitude
          0.801           2.661           6.236   for wavelenght range in A
minimization took 8.29712605476 seconds
```
The main plot out will look like this:
![alt tag](https://raw.githubusercontent.com/nyusngroup/SESNspectraLib/master/10qts_20100815_Lick_3-m_v1-z.flm-flat-Fe.png)



###Required packages:
This code makes use of several default python packages, including [NumPy](http://www.numpy.org/), [SciPy](https://www.scipy.org/), [Matplotlib](http://matplotlib.org/) and [pickle](https://docs.python.org/2/library/pickle.html), and a few additional, less common, packages:

- [emcee](http://dan.iel.fm/emcee/current/)
- [corner](https://github.com/dfm/corner.py) (optional, but recommanded: it generates plots of the marginalized distributions as visual diagnostics of your fit)

If you use data products in this repository, please <b>cite</b> the following references. Here are the bibtex entries for your convenience, as well as the DOI from Zenodo for this repository (TO BE FINISHED).

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
