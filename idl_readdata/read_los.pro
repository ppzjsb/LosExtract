
;;-------------------------------------------------------------------

;; Read LOS and tau files

;;-------------------------------------------------------------------

ztime  = 0.0d
omega0 = 0.0d
omegaL = 0.0d
omegab = 0.0d
h100   = 0.0d
box100 = 0.0d
Xh     = 0.0d
nbins  = 0L
numlos = 0L

openr,1,filename1

;; Header
readu,1,ztime,omega0,omegaL,omegab,h100,box100,Xh
readu,1,nbins,numlos

;; Sight-line locations (axis, x,y,z)
iaxis = lonarr(numlos)
xaxis = dblarr(numlos)
yaxis = dblarr(numlos)
zaxis = dblarr(numlos)
readu,1,iaxis,xaxis,yaxis,zaxis


;; Sight-line scale (comoving kpc/h, km/s)
posaxis = dblarr(nbins)
velaxis = dblarr(nbins)
readu,1,posaxis,velaxis

;; Normalised gas density, rho/<rho>
density = dblarr(numlos*nbins)  
readu,1,density

;; HI data (HI/H, T in K, vpec in km/s)
H1frac   = dblarr(numlos*nbins) 
temp_H1  = dblarr(numlos*nbins)
vpec_H1  = dblarr(numlos*nbins)
readu,1,H1frac,temp_H1,vpec_H1

if HYBRID_RT eq 1 then begin
   gJH0 = dblarr(numlos*nbins)
   readu,1,gJH0
endif

if HE2_FLAG eq 1 then begin
;; HeII data (HeIII/H, T in K, vpec in km/s)
   He2frac   = dblarr(numlos*nbins) 
   temp_He2  = dblarr(numlos*nbins)
   vpec_He2  = dblarr(numlos*nbins)
   readu,1,He2frac,temp_He2,vpec_He2
endif

close,1


;; HI Lya optical depth
openr,2,filename2
tau_H1 = dblarr(numlos*nbins)
readu,2,tau_H1
close,2

;; Check for NaNs
ind = where(finite(tau_H1) eq 0)
if(ind(0) ne -1) then begin
   print,'Bad elements in tau_H1!',n_elements(ind)
   print,ind(0)/nbins
endif


;; Optionally read optical depth quantities 
if TAUW_FLAG eq 1 then begin
   openr,3,filename3
   density_tau  = dblarr(numlos*nbins) 
   temp_H1_tau  = dblarr(numlos*nbins)
   readu,3,density_tau,temp_H1_tau
   close,3   
endif

if HE2_FLAG eq 1 then begin
   
;; HeII Lya optical depth
   openr,4,filename4
   tau_He2 = dblarr(numlos*nbins)
   readu,4,tau_He2
   close,4
   
;; Check for NaNs
   ind = where(finite(tau_He2) eq 0)
   if(ind(0) ne -1) then begin
      print,'Bad elements in tau_He2!',n_elements(ind)
   endif
   
;; Optionally read optical depth quantities 
   if TAUW_FLAG eq 1 then begin
      openr,5,filename5
      density_tau_He2 = dblarr(numlos*nbins) 
      temp_He2_tau    = dblarr(numlos*nbins)
      readu,5,density_tau_He2,temp_He2_tau
      close,5   
   endif

endif


if SILICON eq 1 then begin
   openr,6,filename6
   tau_Si2_1190 = dblarr(numlos*nbins)
   readu,6,tau_Si2_1190
   close,6
   
   openr,7,filename7
   tau_Si2_1193 = dblarr(numlos*nbins)
   readu,7,tau_Si2_1193
   close,7
   
   openr,8,filename8
   tau_Si2_1260 = dblarr(numlos*nbins)
   readu,8,tau_Si2_1260
   close,8
   
   openr,9,filename9
   tau_Si3_1207 = dblarr(numlos*nbins)
   readu,9,tau_Si3_1207
   close,9

   openr,10,filename10
   tau_Si4_1394 = dblarr(numlos*nbins)
   readu,10,tau_Si4_1394
   close,10

   openr,11,filename11
   tau_Si4_1403 = dblarr(numlos*nbins)
   readu,11,tau_Si4_1403
   close,11
   
endif


;;-------------------------------------------------------------------
