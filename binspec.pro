;;
;; binspec - simple binning of spectra 
;;

function binspec, wavelength, flux, wstart, wend, wbin, outlam

;
; Simple binning of spectra
; (sent to me by Gregory Sainton)
;  
; Syntax - outflux = binspec(wavelength, flux, wstart, wend, wbin, outlam)
;
; wavelength [%f] - input wavelength vector
; flux       [%f] - corresponding input flux vector
; wstart     [%f] - starting wavelength of output wavelength vector
; wend       [%f] - ending wavelength of output wavelength vector
; wbin       [%f] - bin size in units of wavelength
; outlam     [%f] - output wavelength vector
;
; Dependencies: function integ.pro [MERON library]
;

nlam   = (wend-wstart)/wbin+1
outlam = findgen(nlam)*wbin+wstart
answer = fltarr(nlam)
 
interplam = [wavelength,outlam]
interplam = interplam(sort(interplam))
interplam = interplam(uniq(interplam))
 
interpflux = interpol(flux,wavelength,interplam)
 
for i=0l,nlam-2 do begin
       w = where(interplam ge outlam[i] and interplam le outlam[i+1])
       if n_elements(w) eq 2 then $
           answer[i]=0.5*(total(interpflux[w])*wbin) else $
           answer[i]=integ(interplam[w],interpflux[w],/val)
endfor
 
answer[nlam-1] = answer[nlam-2]

answer[where(outlam GE max(wavelength) or outlam LT min(wavelength))]=0  ;add by Yuqian

return, answer/wbin

end




