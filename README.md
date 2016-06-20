# SESNspectraLib

This repository contains code that was used in [Liu et al. (2016)](http://arxiv.org/abs/1510.08049) and [Modjaz et al. (2016)](http://arxiv.org/abs/1509.07124) to process Stripped Envelope Supenova (SESNe) spectra. 
Common issues encountered in using SN spectral data include missing uncertainties and inhomogeneous measurements of the velocity for differen SN types. Our library inteds to address these issues. 

At this time two pieces of code are provided: 
- **FFT_smooth.pro**: IDL routine to smooth SN spectra by separating SN signal from noise in Fourier space. It should be used to evaluate the uncertainty of spectra in absence of reduction-generated spectral uncertainty arrays, or to reaplce those in a comprehensive manner for a large and diverse dataset [Liu et al. (2016)](http://arxiv.org/abs/1510.08049).
- **Ic_conv_Icbl_MCMC.py**: Python module to measure the absorption velocity of SN Ic-bl from the Fe I 5169 spectral feature. This code provides and implements a method to measure the spectral velocity of SN Ic-bl consistently with what is usually done for oher subtypes, for fair and meanungful comparisons, enabling the measurement of the Fe I 5169 blue-shift despite the fact that this feature is blended for Ic-bl. The code measures the blue-shift of Fe I 5169 blue-shift by matching the SN Ic-bl spectrum under consideration to a broadened and blueshifted SN Ic template, at similar phases. We use it for Fe I 5169 and [Modjaz et al. (2016)](http://arxiv.org/abs/1509.07124) but it can be used for other features as well.


NOTE: When not specified otherwise (e.g. Ic_conv_Icbl) our code can be used on any medium resolution spectra of SNe of any type, although it is deisgned for and tested on SESN.


If you use data products in this repository, please <b>cite</b> related references listed above.
