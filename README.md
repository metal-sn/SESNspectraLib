# SESNspectraLib

This repository contains code that was used in [Liu et al. (2016)](http://arxiv.org/abs/1510.08049) and [Modjaz et al. (2016)](http://arxiv.org/abs/1509.07124) to process Stripped Envelope Supenova (SESNe) spectra. 
Common issues encountered in using SN spectral data include missing uncertainties and inhomogeneous measurements of the velocity for differen SN types. Our library intends to address these issues. 

At this time two pieces of code are provided: 
- **FFT_smooth.pro**: IDL routine to smooth SN spectra by separating SN signal from noise in Fourier space. It should be used to evaluate the uncertainty of spectra in absence of reduction-generated spectral uncertainty arrays, or to reaplce those in a comprehensive manner for a large and diverse dataset ([Liu et al. (2016)](http://arxiv.org/abs/1510.08049)).
- **Ic_conv_Icbl_MCMC.py**: Python module to measure the absorption velocity of SN Ic-bl from the Fe I 5169 spectral feature. This code provides and implements a method to measure the spectral velocity of SN Ic-bl consistently with what is usually done for oher subtypes, for fair and meanungful comparisons, enabling the measurement of the Fe I 5169 blue-shift despite the fact that this feature is blended in Ic-bl specra (typically the Fe I triplet appears as a single blended feature in Ic-bl). The code measures the blue-shift and width of Fe I 5169 by matching the SN Ic-bl spectrum under consideration to a broadened and blueshifted SN Ic mean template at similar phases. We generally use it for Fe I 5169 ([Modjaz et al. (2016)](http://arxiv.org/abs/1509.07124)) but it can be used for other features as well.


NOTE: When not specified otherwise (e.g. Ic_conv_Icbl) our code can be used on any medium resolution spectra of SNe of any type, although it is deisgned for and tested on SESN.

### Examples
Following code is an example of **FFT_smooth.pro** in IDL, which will produce the top panel of Figure 17 in [Liu et al. (2016)](http://arxiv.org/abs/1510.08049):
```
IDL> readcol, 'sn2004gq-20041212-z-bl.flm', w, f
IDL> FFT_smooth, w, f, 1000, w_ft, f_ft, sep_vel
IDL> plot, w, f
IDL> oplot, w_ft, f_ft, color=250
```
![alt tag](https://github.com/nyusngroup/SESNspectraLib/blob/master/example_IDL_plot.png)

Following code is an example of **Ic_conv_Icbl_MCMC.py** in Python, which will produce: 
- **10qts_20100815_Lick_3-m_v1-z.flm-flat.sav-Fe.dat**, a text file that contains initial template fit region, mean acceptance fraction, initial values for parameters, and 16th, 50th, 84th percentiles of marginalized distribution of model parameters.
- **10qts_20100815_Lick_3-m_v1-z.flm-flat.sav-Fe.p**, a pickle file that contains marginalized distribution of model parameters.
- **10qts_20100815_Lick_3-m_v1-z.flm-flat.sav-Fe.pdf**, a PDF file that contains template fit plot, corner plot, and chain plot.
```
>>> import Ic_conv_Icbl_MCMC as Ic_conv_Icbl
>>> Ic_conv_Icbl.conv('10qts_20100815_Lick_3-m_v1-z.flm-flat.sav', 'meanspecIc_0.sav')
```


If you use data products in this repository, please <b>cite</b> related references listed above.
