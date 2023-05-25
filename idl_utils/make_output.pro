
;;-------------------------------------------------------------------

;; Simple routine to make outputs.txt list for Gadget-3/4

;;-------------------------------------------------------------------

pro make_output

;; Flag=0, generate from list.  Flag != 0, generate equally spaced
;; outputs in time
Flag = 0

filename = "outputs.txt"
openw,1,filename



if (Flag eq 0) then begin

;; Flag = 0.  Modify here
;;-------------------------------------------------------------------
   
   ;; Sherwood outputs
   z = [10.0d, 8.0d, 7.0d, 6.0d, 5.4d, 4.8d, 4.2d, 3.6d, 3.2d, 2.8d, 2.4d, 2.0d]
   
;;-------------------------------------------------------------------

   
   print
   print,'Number of snapshots:',n_elements(z)
   print
   atime = 1.0d/(1.0d + z)
   

   for i=0, n_elements(z) - 1 do begin
      print,'[z, a(t)]',z[i],atime[i]
      printf,1,atime(i)
   endfor
   print

   
endif else begin
   
   
;; Flag !=0.  Modify here
;;-------------------------------------------------------------------
   
   ZSTART = 20.0d
   ZEND   = 4.0d
   DT_MYR = 40.0  ;; Myr
   
   OmegaM = 0.308d
   OmegaL = 0.692d
   h100   = 0.678d

;;-------------------------------------------------------------------
   
   
   MPC = 3.085677581d24
   H   = 1.0d7 / MPC  
   dt  = DT_MYR * 1.0d6 * 3600.0 * 24.0 * 365.0
   
   ;; t(z), e.g. Padmanabhan, Theoretical Astrophysics Vol 3., p185,
   ;; Eq. 3.120
   atime   = 1.0d / (1.0d + ZSTART)
   t_start = 2.0d / (3.0d * H * h100 * sqrt(OMEGAL)) * asinh( sqrt(OMEGAL/OMEGAM) * atime^1.5d )
   
   atime   = 1.0d / (1.0d + ZEND)   
   t_end   = 2.0d / (3.0d * H * h100 * sqrt(OMEGAL)) * asinh( sqrt(OMEGAL/OMEGAM) * atime^1.5d ) 
   
   N_outputs = ceil((t_end - t_start) / dt)
   print
   print,'Number of snapshots:',N_outputs
   print
   
   for i=0, N_outputs do begin
      atime  = ( sqrt(OMEGAM/OMEGAL) * sinh(1.5d * H * h100 * sqrt(OMEGAL) * t_start) )^(2.0d/3.0d) 
      z      = 1.0/atime - 1.0                                                     
      print,'[z, a(t)]',z,atime
      printf,1,atime
      t_start += dt
   endfor
   print

endelse

close,1

end
    
;;-------------------------------------------------------------------
