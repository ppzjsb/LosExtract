
;;-------------------------------------------------------------------

;; Read in and plot one of the pre-computed UV background files

;;-------------------------------------------------------------------

pro plot_UVB

base = './'

filename = base+'./TREECOOL_P19'
openr,1,filename

redshift = dblarr(1000)
gH0      = dblarr(1000)   
eH0      = dblarr(1000)   
gHe0     = dblarr(1000)   
eHe0     = dblarr(1000)
gHep     = dblarr(1000)   
eHep     = dblarr(1000)

;; Read in data
count = 0   
while (not EOF(1)) do begin
   readf,1,in1,in2,in3,in4,in5,in6,in7
   redshift(count) = 10.0^in1 - 1.0d
   gH0(count)  = in2
   eH0(count)  = in5 
   gHe0(count) = in3
   eHe0(count) = in6
   gHep(count) = in4
   eHep(count) = in7
   
   count = count+1
endwhile
close,1

redshift = redshift[0:count-1]
gH0  = gH0[0:count-1]
eH0  = eH0[0:count-1] 
gHe0 = gHe0[0:count-1]
eHe0 = eHe0[0:count-1]
gHep = gHep[0:count-1]
eHep = eHep[0:count-1]


print 
print,filename
print,'UVB file entries:',n_elements(redshift)
print

;;-------------------------------------------------------------------


window,0,xsize=500,ysize=500
Device,Retain=2,true_color=24,decomposed=0
!p.font=-1
plot, redshift, alog10(gH0),linestyle=0,xstyle=1,ystyle=1,xrange=[0,16],yrange=[-21.0,-11.5],xtitle ='z',ytitle='log(photoionisation rate/s^-1)',charsize=1.75,/noerase
loadct,5
oplot,redshift,alog10(gHe0),linestyle=0,color=50
oplot,redshift,alog10(gHep),linestyle=2,color=100

eV = 1.602176565e-12

window,1,xsize=500,ysize=500
Device,Retain=2,true_color=24,decomposed=0
!p.font=-1
plot, redshift, alog10(eH0/gH0/eV),linestyle=0,xstyle=1,ystyle=1,xrange=[0,16],yrange=[0,3],xtitle ='z',ytitle='log(Energy per ionisation/eV)',charsize=1.75,/noerase
loadct,5
oplot,redshift,alog10(eHe0/gHe0/eV),linestyle=0,color=50
oplot,redshift,alog10(eHep/gHep/eV),linestyle=2,color=100

;;-------------------------------------------------------------------


end
