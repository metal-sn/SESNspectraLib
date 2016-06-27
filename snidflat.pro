function snidflat, w, f, z, $ ;input wavelength and flux arrays; input redshift
                   wlog=wlog, fbin=fbin, fmean=fmean, ffilt=ffilt, logspl=logspl, $ ;optional outputs
                   wmin=wmin, wmax=wmax, $ ;restrict *observed* wavelength range of input
                   w0=w0, w1=w1, nw=nw, $ ;[w0,w1] and number of wavelength bins for log wavelength array
                   filter=filter, $ ;set /filter for filtering
                   k1=k1, k2=k2, k3=k3, k4=k4, $ ;wavenumbers for filtering
                   apod=apod, $ ;set apod to percent to apodize
                   fnu=fnu, flambda=flambda, $ ;snidbin options
                   exit=exit      ;returns >0 integer if successful, <0 otherwise

;
; Flatten input spectrum as in SNID
;

exit = -1
fflat = f

;default parameters
if not keyword_set(wmin) then wmin = min(w)
if not keyword_set(wmax) then wmax = max(w)
if not keyword_set(w0)   then w0 = 2500.
if not keyword_set(w1)   then w1 = 10000.
if not keyword_set(nw)   then nw = 1024
if not keyword_set(k1)   then k1 = 1
if not keyword_set(k2)   then k2 = 4
if not keyword_set(k3)   then k3 = nw / 12
if not keyword_set(k4)   then k4 = nw / 10
if not keyword_set(apod) then apod = 5.0

;set up log wavelength scale
dwlog = alog(w1/w0)/nw
;wlog = fltarr(nw+1)
;for i=0, nw do wlog[i]=w0*exp(float(i)*dwlog)
wlog = w0*exp(indgen(nw+1)*dwlog)
for i=0, nw-1 do wlog[i]=0.5*(wlog[i]+wlog[i+1]) ;half-wavelength bins (cf. snidbin)
;wlog = wlog[0:nw-1]

;de-redshift input spectrum
wz = w/(1.+z)

;rebin input onto log wavelength scale
rr = where(w ge wmin and w le wmax,nrr) ;wmin,wmax apply to *observed* wavelength range!
fbin = snidbin(wz[rr],f[rr],nw,w0,dwlog,fnu=keyword_set(fnu),flambda=keyword_set(flambda)) ;bin using *de-redshifted* wavelength array!

;remove the mean with a spline fit
izoff = -1
fmean = meanzero(l1,l2,fbin,izoff,nknod,xknod,yknod,logspl=logspl)
if (n_elements(fmean) eq 1 and fmean[0] eq -1) then return, -1

;apply bandpass filter
ffilt = fmean
if keyword_set(filter) then begin
    ftmp=snidfilt(fmean,k1=k1,k2=k2,k3=k3,k4=k4)    
    ffilt[l1:l2]=ftmp[l1:l2]
endif

;apodize the end
fflat = ffilt
if keyword_set(apod) then $
  fflat=apodize(l1,l2,ffilt,apod)

exit = 1
return, fflat
end
