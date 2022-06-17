
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
close,1

;; HI Lya optical depth
openr,2,filename2
tau_H1 = dblarr(numlos*nbins)
readu,2,tau_H1
close,2

flux = dblarr(numlos*nbins)  

;; Check for NaNs
ind = where(finite(tau_H1) eq 0)
if(ind(0) ne -1) then begin
   print,'Bad elements in tau_H1!',n_elements(ind)
endif


;; Optionally read optical depth quantities 
if TAUW_FLAG eq 1 then begin
   openr,3,filename3
   density_tau  = dblarr(numlos*nbins) 
   temp_H1_tau  = dblarr(numlos*nbins)
   readu,3,density_tau,temp_H1_tau
   close,3   
endif

;;-------------------------------------------------------------------
