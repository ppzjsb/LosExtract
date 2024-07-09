
;;-------------------------------------------------------------------

;; Read in the output of LosExtract and plot some quantites.

;;-------------------------------------------------------------------

pro get_losdata

;; Flag for optical depth weighted quantities (0 or 1)
TAUW_FLAG = 0

;; Flag for HeII Lyman-alpha (0 or 1)
HE2_FLAG  = 0

;; Flag for SiII and SiIII absorption (0 or 1)
SILICON   = 1

;; Flag for test of line profile convolution (0 or 1).  The GP optical
;; depths should be recovered.
TEST_KERNEL = 0

;; Select LOS (0 to NUMLOS-1)
PLOTLOS = 42

base      = '../../'

;;-------------------------------------------------------------------

;;filename1 = base+'spec1024_n5000_z2.000.dat'
filename1 = base+'los2048_n5000_z3.000.dat'
filename2 = base+'tauH1_v2048_n5000_z3.000.dat'

if TAUW_FLAG eq 1 then  begin
   filename3 = base+'tauwH1_x256_n5000_z2.000.dat'
endif

if HE2_FLAG eq 1 then begin
   filename4 = base+'tauHe2_x256_n5000_z2.000.dat'
   
   if TAUW_FLAG eq 1 then  begin
      filename5 = base+'tauwHe2_x256_n5000_z2.000.dat'
   endif
endif

if SILICON eq 1 then begin
   filename6 = base+'tauSi2_1190_v2048_n5000_z3.000.dat'
   filename7 = base+'tauSi2_1193_v2048_n5000_z3.000.dat'
   filename8 = base+'tauSi2_1260_v2048_n5000_z3.000.dat'
   filename9 = base+'tauSi3_1207_v2048_n5000_z3.000.dat'
endif
   
@read_los
;;@read_spec

@constants

dv = velaxis[1]-velaxis[0]

ind = where(density gt 1000.0 and temp_H1 lt 1.0e5)
if(ind(0) ne -1) then begin
   print
   print,'Pixels where Delta > 1000 and log*T/K) < 5:',n_elements(ind)
   print,alog10(temp_H1(ind))
   print,alog10(density(ind))
endif
   

print
print,filename1
print,filename2
print,'z:        ',ztime
print,'nbins:    ',nbins
print,'numlos:   ',numlos
print,'OmegaM:   ',omega0
print,'OmegaL:   ',omegaL
print,'Omegab:   ',omegab
print,'Hubble:   ',h100
print,'BoxSize:  ',box100
print,'Xh:       ',Xh
print,'<F>:      ',mean(exp(-tau_H1[*]))
print,'dv [km/s]:',dv
print


print,'Range of log(Delta)'
print,alog10(min(density)),alog10(max(density))
print

print,'Range of log(T/K)'
print,alog10(min(temp_H1)),alog10(max(temp_H1))
print

print,'Range of log(tau_H1)'
print,alog10(min(tau_H1)),alog10(max(tau_H1))
print

if HE2_FLAG eq 1 then begin
   print,'Range of log(tau_He2)'
   print,alog10(min(tau_He2)),alog10(max(tau_He2))
   print
endif


;; Plot transmission in first sight-line
window,0,xsize=1200,ysize=350,title='HI-Lya'
Device,Retain=2,true_color=24,decomposed=0
!p.font=-1
plot,velaxis,exp(-tau_H1[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]),xstyle=1,ystyle=1,charsize=1.75,yrange=[-0.1,1.1],ytitle='Transmitted flux (HI-Lya)',xtitle='Hubble velocity [km/s]'
if (TAUW_FLAG eq 1) then begin
   loadct,5
   oplot,velaxis,(alog10(density_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]) - min(alog10(density_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1])))/max(alog10(density_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1])- min(alog10(density_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]))),linestyle=1,color=50
   if(TEST_KERNEL ne 1) then begin
      oplot,velaxis,(alog10(temp_H1_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]) - min(alog10(temp_H1_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1])))/max(alog10(temp_H1_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1])- min(alog10(temp_H1_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]))),linestyle=2,color=100
   endif
endif


if HE2_FLAG eq 1 then begin
   window,1,xsize=1200,ysize=350,title='HeII-Lya'
   Device,Retain=2,true_color=24,decomposed=0
   !p.font=-1
   plot,velaxis,exp(-tau_He2[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]),xstyle=1,ystyle=1,charsize=1.75,yrange=[-0.1,1.1],ytitle='Transmitted flux (HeII-Lya)',xtitle='Hubble velocity [km/s]'
   if (TAUW_FLAG eq 1) then begin
      loadct,5
      oplot,velaxis,(alog10(density_tau_He2[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]) - min(alog10(density_tau_He2[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1])))/max(alog10(density_tau_He2[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1])- min(alog10(density_tau_He2[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]))),linestyle=1,color=50
      if(TEST_KERNEL ne 1) then begin
         oplot,velaxis,(alog10(temp_He2_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]) - min(alog10(temp_He2_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1])))/max(alog10(temp_He2_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1])- min(alog10(temp_He2_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]))),linestyle=2,color=100
      endif
   endif
endif

if SILICON eq 1 then begin
   window,1,xsize=1200,ysize=350,title='SiII and SiIII'
   Device,Retain=2,true_color=24,decomposed=0
   !p.font=-1
   plot,velaxis,exp(-tau_Si2_1190[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]),yrange=[-0.1,1.1],xstyle=1,ystyle=1,charsize=1.75,ytitle='Transmitted flux (SiII, SiIII)',xtitle='Hubble velocity [km/s]'
   oplot,velaxis,exp(-tau_Si2_1193[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]),linestyle=1,color=150
   oplot,velaxis,exp(-tau_Si2_1260[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]),linestyle=2,color=175
   oplot,velaxis,exp(-tau_Si3_1207[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]),linestyle=3,color=200
endif

;; Plot temperature-density plane, 10000 random points
ind = floor(n_elements(density)*randomu(10, 10000,/uniform))
window,2,xsize=500,ysize=500
Device,Retain=2
!p.font=-1
plot,alog10(density(ind)),alog10(temp_H1(ind)),xstyle=1,ystyle=1,charsize=1.75,xtitle='log Overdensity',ytitle='log Temperature',psym=3,xrange=[-2,3],yrange=[2,6]
if HE2_FLAG eq 1 then begin
   oplot,alog10(density(ind)),alog10(temp_He2(ind)),psym=3,color=100
endif



;;-------------------------------------------------------------------

if(TEST_KERNEL eq 1) then begin
   
;; Optical depth PDF
   binmin  = -8.0
   binmax  = 8.0
   binsize = 0.25
   pdfbins = (binmax-binmin)/binsize 

   ind = where(tau_H1 gt 0.0)
   tau_pdf = histogram(alog10(tau_H1(ind)),MIN=binmin,BINSIZE=binsize,NBINS=pdfbins,locations=logtaubin)   
   tau_pdf /= (n_elements(tau_H1(ind))*binsize)
   
   if HE2_FLAG eq 1 then begin
      ind = where(tau_He2 gt 0.0)
      tau_pdf_He2 = histogram(alog10(tau_He2(ind)),MIN=binmin,BINSIZE=binsize,NBINS=pdfbins,locations=logtaubin)
      tau_pdf_He2 /= (n_elements(tau_He2(ind))*binsize)
   endif

   logtaubin += 0.5*binsize
   
   loadct,5
   window,3,xsize=500,ysize=500
   Device,Retain=2,true_color=24,decomposed=0
   !p.font=-1   
   plot,logtaubin,tau_pdf,xstyle=1,ystyle=1,charsize=1.75,/yl,xrange=[binmin,binmax],yrange=[1.0d-7,10],xtitle='p(log Optical depth)',ytitle='log(Optical depth)'
   oplot,logtaubin,tau_pdf,psym=4,color=140,symsize=1.2
   if HE2_FLAG eq 1 then begin
      oplot,logtaubin,tau_pdf_He2,linestyle=1
      oplot,logtaubin,tau_pdf_He2,psym=4,symsize=1.2
   endif
   
   
;; Contour plot of optical depth against log(Delta)
   binx = 0.05d
   minx = -3.0d - binx*0.5d
   maxx = 4.0d  + binx*0.5d
   
   biny = 0.05d
   miny = -8.0d - biny*0.5d
   maxy = 8.0d + biny*0.5d
   
   xrow  = floor((maxx - minx)/binx)+1
   yrow  = floor((maxy - miny)/biny)+1
   xcont = dblarr(xrow,yrow)
   ycont = dblarr(xrow,yrow)
   
   for j=0,yrow-1 do begin
      for i=0,xrow-1 do begin 
         kk = i + xrow*j  
         xcont[i,j] = minx + i*binx + binx*0.5d
         ycont[i,j] = miny + j*biny + biny*0.5d
      endfor
   endfor

   zcont = HIST_2D(alog10(density),alog10(tau_H1),bin1=binx,bin2=biny,max1=maxx,max2=maxy,min1=minx,min2=miny)

   if HE2_FLAG eq 1 then begin
      zcont_He2 = HIST_2D(alog10(density),alog10(tau_He2),bin1=binx,bin2=biny,max1=maxx,max2=maxy,min1=minx,min2=miny)
   endif
  
   
;; Compare to Gunn-Peterson optical depth for a consistency test.
;; Should match exactly for an isothermal IGM with constant H1 or He2
;; fraction and ignoring peculiar velocities.  If it does not match,
;; the broadening kernel is under-resolved and the convolution needs
;; to be performed on a finer grid.  Note that LosExtract has an
;; option, RESAMPLE_KERNEL that should ensure it works.
   
   
   log_Delta = -4.0d + dindgen(1000)/100.0
   hubble    = 1.0d7 / MPC
   yhelium   = (1.0d - Xh)/(4.0d * Xh)       
   
   rhoc   = 3.0d * (h100*hubble)^2.0/(8.0d * !dpi * GRAVITY) 
   rho    = rhoc * omegab * (1.0d + ztime)^3.0d
   nHcgs  = rho * XH / PROTONMASS
   Hz     = (1.0d7/MPC) * h100 * sqrt(omega0 * (1.0d + ztime)^3.0d + omegaL)

   xH1 = 1.0e-5 ;; H1/H
   sigma_Lya = sqrt(3.0*!dpi*SIGMA_T/8.0) * LAMBDA_LYA_H1 * FOSC_LYA
   tau_GP = sigma_Lya * C_CMS * xH1 * nHcgs * (10.0^log_Delta) / Hz
   
   window,4,xsize=500,ysize=500
   Device,Retain=2,true_color=24,decomposed=0
   !p.font=-1
   contour,zcont,xcont,ycont,/fill,xtitle='log Overdensity',ytitle='log Optical depth',xstyle=1,ystyle=1,charsize=1.75,xrange=[-2,3],yrange=[-6,6],levels=[10.0^0.0,10.0^0.5,10.0^1,10.0^1.5,10.0^2,10.0^2.5,10.0^3,10.0^3.5,10.0^4,10.0^4.5,10.0^5.0]
   oplot,log_Delta,alog10(tau_GP),linestyle=0,color=140
   
   if HE2_FLAG eq 1 then begin
      xHe2 = 1.0e-4 ;; He2/H
      sigma_Lya_He2 = sqrt(3.0*!dpi*SIGMA_T/8.0) * LAMBDA_LYA_HE2 * FOSC_LYA
      tau_GP_He2 = sigma_Lya_He2 * C_CMS * xHe2 * nHcgs * (10.0^log_Delta) / Hz
      
      loadct,0
      contour,zcont_He2,xcont,ycont,/fill,xtitle='log Overdensity',ytitle='log Optical depth',xstyle=1,ystyle=1,charsize=1.75,xrange=[-2,3],yrange=[-6,6],levels=[10.0^0.0,10.0^0.5,10.0^1,10.0^1.5,10.0^2,10.0^2.5,10.0^3,10.0^3.5,10.0^4,10.0^4.5,10.0^5.0],/overplot
      oplot,log_Delta,alog10(tau_GP_He2),linestyle=1
   endif
endif

;;-------------------------------------------------------------------

end
