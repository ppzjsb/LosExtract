
;;-------------------------------------------------------------------

;; Rescale mean transmission/effective optical depth 
      if (RESCALE eq 1) then begin
         
;; Becker et al. 2013, MNRAS , 430, 2067, Eq. 5 to 2.2<=z < 4.4 ;;
         if(ztime lt 4.4) then begin
            taueff = -0.132d + 0.751d * ((1.0d + ztime)/4.5d)^2.90d 
         endif
         
;; Similar to Viel et al. (2013), PRD, 88, 043502, Eq. 4, 4.4<= z < 5.5 only 
         if(ztime ge 4.4) then begin
            taueff = 1.142d * ((1.0d + ztime)/5.4d)^4.91d
         endif
         
         flux_obs = exp(-taueff)
         scale = 0.01d
         for iscale=1,1000 do begin  ;; This is just Newton-Raphson 
            scale_old = scale
            scale = scale +  ( mean(exp(-scale * tau_lya_H1[*])) - flux_obs ) $
                    / mean(tau_lya_H1[*] * exp(-scale * tau_lya_H1[*]))
            dscale = abs(scale - scale_old)
            if (dscale lt 1.0d-8) then break
         endfor
         flux = exp(-1.0d * scale * tau_lya_H1)
      endif else begin
         flux = exp(-1.0d * tau_lya_H1)
      endelse
      
      ;; Get the box size in km/s
      atime  = 1.0d/(1.0d + ztime)
      Hz     = 100.0d *h100*sqrt(omega0/(atime*atime*atime)+omegaL) ;; /* km s^-1 Mpc^-1 */
      boxkms = Hz * box100/(1.0d3*h100*(1.0d + ztime))

      
      ;; Use a Gaussian with specified FWHM
      if (CONVOLVE eq 1) then begin
         
         FWHM    = 7.0d ;; km/s
         
         FWHMrel = FWHM * nbins/boxkms 
         sigma   = FWHMrel/(2.0d * sqrt(2.0d * alog(2.0d)))
         
         xx     = dindgen(nbins) - nbins/2.0d
         kernel = 1.0d/(sigma * sqrt(2.0d * !dpi))  * exp(-0.5d * (xx/sigma)*(xx/sigma))
         
         for iconv=0,numlos-1 do begin
            flux[iconv*nbins :(iconv+1)*nbins-1] = $
               convol(flux[iconv*nbins :(iconv+1)*nbins-1],kernel,/edge_wrap,/normalize)
         endfor
      endif
      
;;-------------------------------------------------------------------
