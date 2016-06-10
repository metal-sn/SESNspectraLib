pro discreteFT_0703, w, f, noise_cut_vel, w_ln_bin, f_bin, num_bin,f_bin_ft_smooth_inv,num_try_vel

; FFT smoothing 
; w and f are original wavelength and flux
; w_ln_bin and f_bin are wavelength and flux used to FFT
; num_bin is the number in w_ln_bin or f_bin
; f_bin_ft_smooth_inv is FFT smoothed flux
; num_try_vel is noise cutting velocity

  w_ln=alog(w)                  ; convert to log(w) space
  num=n_elements(w_ln)
  binsize_noise=w_ln[1]-w_ln[0] ; largest width of noise      
  binsize = (w_ln[num-1]-w_ln[num-2])
      ;print, 'largest velocity shift (km/s) in ln(w): ',binsize_noise,binsize_noise*3.0e5
      ;print, 'smallest velocity shift (km/s) in ln(w): ',binsize, binsize*3.0e5
      f_bin = binspec(w_ln,f,min(w_ln),max(w_ln),binsize,w_ln_bin) ;bin to same sampling interval in log(w) space
      f_bin=f_bin(where(f_bin NE 0))
      w_ln_bin=w_ln_bin(where(f_bin NE 0))
      num_bin=n_elements(f_bin)
      ;print, 'number of elements before binning: ',num
      ;print, 'number of elements after binning: ', num_bin

; signal to noise
      ;snr[n]=der_snr(f_bin)

; flourier transform
      f_bin_ft=fft(f_bin,-1)*num_bin
      ;print, f_bin_ft(1)
      f_bin_ft_inv=float(fft(f_bin_ft,1)) ; just check if inverse FT works well
;calculate frequency
      X = (FINDGEN((num_bin - 1)/2) + 1)
      is_N_even = (num_bin MOD 2) EQ 0
      if (is_N_even) then $
         freq = [0.0, X, num_bin/2, -num_bin/2 + X]/num_bin $
      else $
         freq = [0.0, X, -(num_bin/2 + 1) + X]/num_bin
      ;print, 'number of fluxes after FT: ', n_elements(f_bin_ft)
      ;print, 'number of frequency in FT: ',n_elements(freq)

;filtering spectra
      
      num_upper=max(where(1.0/freq[0:num_bin-1]*3.0e5*binsize GT noise_cut_vel))
      num_lower=max(where(1.0/freq[0:num_bin-1]*3.0e5*binsize GT 1.0e5, num_num_lower))
      f_bin_ft_line=fltarr(num_upper)
      if num_num_lower lt 3 then num_lower=min([num_bin/100,10])
      intercept=mean(abs(f_bin_ft(num_lower :num_upper)))
      ;print, num_upper
      ;print, intercept
      f_bin_ft_line=(f_bin_ft_line+1.0)*intercept ; a straight line along x axis

      coeffs1=powerlaw_fit(freq(num_lower :num_bin/2), abs(f_bin_ft(num_lower :num_bin/2)))
      
      ;print, 'parameters of a power law fitting: ', coeffs1
      ;print, 'first several amplitude of FT: ',abs(f_bin_ft(0:5))
      delta=(freq(1:num_bin/2)/coeffs1[0])^coeffs1[1]-intercept ;-sigma_use

      if (max(delta) lt 0) and (num_upper*2 le num_bin/2) then begin
; initial guess      
         g = linear_fit(alog(freq[num_lower :num_upper*2]), alog(f_bin_ft[num_lower :num_upper*2]))
         a = [10e1^(-g[0]/g[1]), g[1]]

         coeffs1=powerlaw_fit(freq(num_lower :num_bin/2), abs(f_bin_ft(num_lower :num_bin/2)),guess=a)
         delta=(freq(1:num_bin/2)/coeffs1[0])^coeffs1[1]-intercept 
         print, 'll'
      endif

      if max(delta) lt 0 then begin
; initial guess      
         g = linear_fit(alog(freq[num_lower*2:num_upper]), alog(f_bin_ft[num_lower*2:num_upper]))
         a = [10e1^(-g[0]/g[1]), g[1]]

         coeffs1=powerlaw_fit(freq(num_lower :num_bin/2), abs(f_bin_ft(num_lower :num_bin/2)),guess=a)
        print, 'llll'
         delta=(freq(1:num_bin/2)/coeffs1[0])^coeffs1[1]-intercept 
      endif

      if max(delta) lt 0 then begin
; initial guess      
         g = linear_fit(alog(freq[num_lower :num_upper]), alog(f_bin_ft[num_lower :num_upper]))
         a = [10e1^(-g[0]/g[1]), g[1]]

         coeffs1=powerlaw_fit(freq(num_lower :num_bin/2), abs(f_bin_ft(num_lower :num_bin/2)),guess=a)
         print, 'llllll'
         delta=(freq(1:num_bin/2)/coeffs1[0])^coeffs1[1]-intercept 
      endif

      if max(delta) lt 0 then begin
; initial guess      
         g = linear_fit(alog(freq[num_lower :10*num_lower]), alog(f_bin_ft[num_lower :10*num_lower]))
         a = [10e1^(-g[0]/g[1]), g[1]]

         coeffs1=powerlaw_fit(freq(num_lower :num_bin/2), abs(f_bin_ft(num_lower :num_bin/2)),guess=a)
         print, 'llllllll'
         delta=(freq(1:num_bin/2)/coeffs1[0])^coeffs1[1]-intercept 
      endif

      if max(delta) lt 0 then begin
         cd,'/home/users/yl1260/regenerateSNID/spec_FFT'
         print, 'try another guess value'
         return
      endif
      num_try=min(where(delta LT 0)) 
      num_try_vel=1.0/freq(num_try)*3.0e5*binsize

      f_bin_ft_smooth=f_bin_ft/num_bin
      for j=1L, num_bin do begin
         if j lt num_bin/2+1 then $
            number = j - 1 else $
               number = j - num_bin - 1
         numa = abs(number)
         if numa gt num_try then begin
            f_bin_ft_smooth[j-1]=0
         endif
      endfor     
      f_bin_ft_smooth_inv=float(fft(f_bin_ft_smooth,1))

END
