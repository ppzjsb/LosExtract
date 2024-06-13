
;;-------------------------------------------------------------------

;; Calculate the Lyman-alpha forest power spectrum.  Provides an
;; important consistency test of the simulation pipeline.

;;-------------------------------------------------------------------

pro get_pk
  
;; Flags for processing spectra
RESCALE   = 0  ;; optionally rescale to observed tau_eff
CONVOLVE  = 0  ;; optionally convolve with Gaussian with FWHM=7km/s
REBIN     = 0  ;; optionally rebin/rsample P(k) to be equally spaced log10(k)

;; Location of input files
base = '../../'

;; Select models
;;modval  = ['planck1_40_512_G3/']
modval  = ['']

;; Select redshifts
;;zval    = ['2.000','3.000','4.000','5.000']
zval    = ['3.000']

;; Set range for rebinning here
kmin    = -2.9d
kmax    = -0.3d
binsize = 0.1d
npts    = 1 + abs(kmax-kmin)/binsize
logk_bin = kmin + dindgen(npts)*binsize
k_bin    = 10.0^logk_bin

;;-------------------------------------------------------------------


for m_ind=0, n_elements(modval)-1 do begin
   for z_ind=0, n_elements(zval)-1 do begin
   
      @read_los
    ;;  @read_spec
      @process_spectra
      
      print,'Box size [km/s]',boxkms
      print,'<F>    :',mean(flux)
      print,'tau_eff:',-alog(mean(flux))
      
;; Flux estimator
      df = flux/mean(flux) - 1.0d
      
      ;;  See IDL help on FFT for further details on why we do this
      X = (findgen((nbins-1)/2) + 1)
      
      if((nbins mod 2) eq 0) then begin ;; case for even N
         pk   = dblarr(nbins/2)
         freq = [0.0, X, nbins/2, -nbins/2 + X]/boxkms
         k    = 2.0d *!pi*freq[1:nbins/2]                  ;; wavenumber, s km^-1
         
         if(z_ind eq 0 and m_ind eq 0) then begin
            flux_pk     = dblarr(n_elements(modval),n_elements(zval),nbins/2)
            flux_pk_bin = dblarr(n_elements(modval),n_elements(zval),npts)
         endif
         
      endif else begin ;; case for odd N
         pk   = dblarr(nbins/2-1)
         freq = [0.0, X, -(nbins/2 +1) + X]/boxkms
         k    = 2.0d *!pi*freq[1:nbins/2-1]
         
         if(z_ind eq 0 and m_ind eq 0) then begin
            flux_pk     = dblarr(n_elements(modval),n_elements(zval),nbins/2-1)
            flux_pk_bin = dblarr(n_elements(modval),n_elements(zval),npts)
         endif
      endelse
      
      
;; Power spectrum is the square of the FT.  Compute for each
;; sight-line individually.
      pk[*] = 0.0
      for iproc=0L,numlos-1L do begin
         ft_field = fft(df(iproc*nbins :(iproc+1)*nbins - 1),/double)
         if((nbins mod 2) eq 0) then begin 
            pk += abs(ft_field[1:nbins/2])*abs(ft_field[1:nbins/2])
         endif else begin 
            pk += abs(ft_field[1:nbins/2-1])*abs(ft_field[1:nbins/2-1])
         endelse
      endfor
      
;; Averaged, dimensionless power spectrum (multiply P_1d(k) by k/!pi)
      pk *= k * boxkms / (!dpi * double(numlos))

     ; print
     ; print,'log k :',alog10(k[0]),alog10(k[nbins/64-1]),alog10(k[nbins/8-1]),alog10(k[nbins/2 - 1])
     ; print,'log pk:',alog10(pk[0]),alog10(pk[nbins/64 - 1]),alog10(pk[nbins/8 - 1]),alog10(pk[nbins/2 - 1])
     ; print
      
      
;; Save in an array for plotting 
      flux_pk[m_ind,z_ind,*] = pk 

      ;; Rebin/resample P(k)   
      if(REBIN eq 1) then begin
         for i=0,npts-1 do begin
            bin_lower = logk_bin[i] - 0.5*binsize
            bin_upper = logk_bin[i] + 0.5*binsize
            
            indbin = where((alog10(k) ge bin_lower) and (alog10(k) le bin_upper))
            
            
            if(indbin(0) eq -1) then begin
               print,'No wavenumbers in bin',bin_lower,bin_upper
               print,'Interpolating...'
               flux_pk_bin[m_ind,z_ind,i] = interpol(pk,alog10(k),logk_bin[i])
            endif else begin
               flux_pk_bin[m_ind,z_ind,i] = mean(pk(indbin))
            endelse
         endfor
      endif
      
   endfor
endfor


;;-------------------------------------------------------------------

window,0,xsize=600,ysize=600
device,retain=2

if(REBIN eq 1) then begin
   plot,logk_bin,alog10(flux_pk_bin[0,0,*]),linestyle=0,xrange=[-3.5,0.5],xstyle=1,yrange=[-20,1],ystyle=1,xtitle='log(k/s/km)',ytitle=' log(kP(k)/pi)',charsize=1.75
;   oplot,logk_bin,alog10(flux_pk_bin[0,1,*]),linestyle=1
;   oplot,logk_bin,alog10(flux_pk_bin[0,2,*]),linestyle=2
;   oplot,logk_bin,alog10(flux_pk_bin[0,3,*]),linestyle=3
endif else begin
   plot,alog10(k),alog10(flux_pk[0,0,*]),linestyle=0,xrange=[-3.5,0.5],xstyle=1,yrange=[-20,1],ystyle=1,xtitle='log(k/s/km)',ytitle=' log(kP(k)/pi)',charsize=1.75
 ;  oplot,alog10(k),alog10(flux_pk[0,1,*]),linestyle=1
 ;  oplot,alog10(k),alog10(flux_pk[0,2,*]),linestyle=2
 ;  oplot,alog10(k),alog10(flux_pk[0,3,*]),linestyle=3
endelse

;;-------------------------------------------------------------------

end
