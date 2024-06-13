
;;-------------------------------------------------------------------

;; Read LOS and tau files

;;-------------------------------------------------------------------


filename1 = base+modval[m_ind]+'los2048_n5000_z'+zval[z_ind]+'.dat'
filename2 = base+modval[m_ind]+'tauH1_v2048_n5000_z'+zval[z_ind]+'.dat'
   
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
print

;; HI Lya optical depth
openr,2,filename2
tau_lya_H1 = dblarr(numlos*nbins)
readu,2,tau_lya_H1
close,2

flux = dblarr(numlos*nbins)  

;; Check for NaNs
ind = where(finite(tau_lya_H1) eq 0)
if(ind(0) ne -1) then begin
   print,'Bad elements in tau_H1!',n_elements(ind)
endif

;;-------------------------------------------------------------------
