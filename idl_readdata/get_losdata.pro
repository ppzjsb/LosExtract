
;;-------------------------------------------------------------------

;; Read in the output of LosExtract and plot some quantites.

;;-------------------------------------------------------------------

pro get_losdata

;; flag for optical depth weighted quantities (0 or 1)
TAUW_FLAG = 0

;; flag for testing the accuracy of the line profile convolution (0 or
;; 1).  Fix xH1 and T to constants by hand in LosExtract, and use
;; NO_PECVEL.  The GP optical depths should be recovered.
TEST_KERNEL = 0

;; Select LOS (0 to NUMLOS-1)
PLOTLOS = 325


base      = '../../../testruns/planck1_40_512_G3/'


;;-------------------------------------------------------------------

;;filename1 = base+'spec1024_n5000_z2.000.dat'
filename1 = base+'los2048_n5000_z5.000.dat'
filename2 = base+'tau2048_n5000_z5.000.dat'

if TAUW_FLAG eq 1 then  begin
   filename3 = base+'tauw1024_n5000_z4.000.dat'
endif

@read_los
;;@read_spec


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
print,'z:      ',ztime
print,'nbins:  ',nbins
print,'numlos: ',numlos
print,'OmegaM: ',omega0
print,'OmegaL: ',omegaL
print,'Omegab: ',omegab
print,'Hubble: ',h100
print,'BoxSize:',box100
print,'Xh:     ',Xh
print,'<F>:    ',mean(exp(-tau_H1[*]))
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


;; Plot transmission in first sight-line
window,0,xsize=1200,ysize=350
 Device,Retain=2,true_color=24,decomposed=0
!p.font=-1
plot,velaxis,exp(-tau_H1[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]),xstyle=1,ystyle=1,charsize=1.75,yrange=[-0.1,1.1],ytitle='Transmitted flux',xtitle='Hubble velocity [km/s]'
if (TAUW_FLAG eq 1) then begin
   loadct,5
   oplot,velaxis,(alog10(density_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]) - min(alog10(density_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1])))/max(alog10(density_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1])- min(alog10(density_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]))),linestyle=1,color=50
   oplot,velaxis,(alog10(temp_H1_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]) - min(alog10(temp_H1_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1])))/max(alog10(temp_H1_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1])- min(alog10(temp_H1_tau[PLOTLOS*nbins :(PLOTLOS+1)*nbins-1]))),linestyle=2,color=100

   
endif


;; Plot temperature-density plane, 10000 random points
ind = floor(n_elements(density)*randomu(10, 10000,/uniform))
window,1,xsize=500,ysize=500
Device,Retain=2
!p.font=-1
plot,alog10(density(ind)),alog10(temp_H1(ind)),xstyle=1,ystyle=1,charsize=1.75,xtitle='log Overdensity',ytitle='log Temperature',psym=3,xrange=[-2,3],yrange=[2,6]



;;-------------------------------------------------------------------

if(TEST_KERNEL eq 1) then begin
   
;; Optical depth PDF
   binmin  = -8.0
   binmax  = 8.0
   binsize = 0.25
   pdfbins = (binmax-binmin)/binsize 

   ind     = where(tau_H1 gt 0.0)
   tau_pdf = histogram(alog10(tau_H1(ind)),MIN=binmin,BINSIZE=binsize,NBINS=pdfbins,locations=logtaubin)
   
   logtaubin += 0.5*binsize
   tau_pdf   /= (n_elements(tau_H1(ind))*binsize)
   
   loadct,5
   window,2,xsize=500,ysize=500
   Device,Retain=2,true_color=24,decomposed=0
   !p.font=-1
   
   plot,logtaubin,tau_pdf,xstyle=1,ystyle=1,charsize=1.75,/yl,xrange=[-10,10],yrange=[1.0d-7,10],xtitle='p(log Optical depth)',ytitle='log(Optical depth)'
   oplot,logtaubin,tau_pdf,psym=4,color=140,symsize=1.2
   

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


;; Compare to Gunn-Peterson optical depth for a consistency test.
;; Should match exactly for an isothermal IGM with constant H1
;; fraction and ignoring peculiar velocities.
   xH1           = 1.0e-5
   GRAVITY       = 6.67428d-8
   MPC           = 3.08568025d24
   PROTONMASS    = 1.672621637d-24
   LAMBDA_LYA_H1 = 1215.6701e-8
   GAMMA_LYA_H1  = 6.265e8 
   
   log_Delta = -4.0d + dindgen(1000)/100.0
   hubble    = 1.0d7 / MPC
   yhelium   = (1.0d - Xh)/(4.0d * Xh)       
   
   rhoc   = 3.0d * (h100*hubble)^2.0/(8.0d * !dpi * GRAVITY) 
   rho    = rhoc * omegab * (1.0d + ztime)^3.0d
   nHcgs  = rho * XH / PROTONMASS
   
   Hz     = (1.0d7/MPC) * h100 * sqrt(omega0 * (1.0d + ztime)^3.0d + omegaL)
   k_Lya  = 3.0d * LAMBDA_LYA_H1^3.0 * GAMMA_LYA_H1/ (8.0d * !dpi * Hz)
   tau_GP = k_Lya * xH1 * nHcgs * (10.0^log_Delta)
   
   
   window,3,xsize=500,ysize=500
   Device,Retain=2,true_color=24,decomposed=0
   !p.font=-1
   contour,zcont,xcont,ycont,/fill,xtitle='log Overdensity',ytitle='log Optical depth',xstyle=1,ystyle=1,charsize=1.75,xrange=[-2.5,3.5],yrange=[-6,6],levels=[10.0^0.0,10.0^0.5,10.0^1,10.0^1.5,10.0^2,10.0^2.5,10.0^3,10.0^3.5,10.0^4,10.0^4.5,10.0^5.0]
   oplot,log_Delta,alog10(tau_GP),linestyle=0,color=100,thick=1

endif

;;-------------------------------------------------------------------

end
