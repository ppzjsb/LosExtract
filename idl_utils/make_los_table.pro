

;;-------------------------------------------------------------------

;; Example of how to make a line-of-sight coordinate table for the
;; SpecExtract LOS_TABLE function

;;-------------------------------------------------------------------

pro make_los_table

;; This outputs the LOS positions used for the Sherwood project.
;; Alternatively, random sight-lines are used.  Note that both these
;; set-ups are already hard coded into SpecExtract. However, more
;; generally this gives one an idea of how to use this function if
;; e.g. extracting through haloes in the FoF catalogue.
  
FLAG_PRACE = 1 

f='table_los_grid.txt'

boxsize = 40.0d ;; Mpc/h


if FLAG_PRACE eq 1 then begin

   ngrid1 = 50L                     
   ngrid2 = 40L                     
   ngrid3 = 30L                     
   
   numlos = ngrid1*ngrid1 + ngrid2*ngrid2 + ngrid3*ngrid3
   iaxis  = intarr(numlos)
   xpos   = dblarr(numlos)
   ypos   = dblarr(numlos)
   zpos   = dblarr(numlos)
   
   
   dgrid1 = boxsize/double(ngrid1)         
   dgrid2 = boxsize/double(ngrid2)
   dgrid3 = boxsize/double(ngrid3)         
   
;; project along x-axis 
   for row=0,ngrid1-1 do begin
      for col=0,ngrid1-1 do begin
         iproc        = row*ngrid1 + col  
         iaxis[iproc] = 1         
         xpos[iproc]  = 0.5d * boxsize
         ypos[iproc]  = (col + 0.131213d)*dgrid1
         zpos[iproc]  = (row + 0.131213d)*dgrid1 
      endfor
   endfor
   
;; project along y-axis 
   for row=0,ngrid2-1 do begin
      for col=0,ngrid2-1 do begin
         iproc        = ngrid1*ngrid1 + row*ngrid2 + col
         iaxis[iproc] = 2         
         ypos[iproc]  = 0.5d * boxsize
         xpos[iproc]  = (col + 0.241008d)*dgrid2
         zpos[iproc]  = (row + 0.241008d)*dgrid2 
      endfor
   endfor
   
;; project along z-axis 
   for row=0,ngrid3-1 do begin
      for col=0,ngrid3-1 do begin
         iproc        = ngrid1*ngrid1 + ngrid2*ngrid2 + row*ngrid3 + col 
         iaxis[iproc] = 3                                     
         zpos[iproc]  = 0.5d * boxsize
         xpos[iproc]  = (col + 0.170482d)*dgrid3        
         ypos[iproc]  = (row + 0.170482d)*dgrid3      
      endfor
   endfor
   
endif else begin
   
   numlos  = 5000L  ;; number of sight-lines
   
   rand1 = randomu(1704, numlos, /uniform,/double)
   rand2 = randomu(1982, numlos, /uniform,/double)
   rand3 = randomu(2410, numlos, /uniform,/double)
   rand4 = randomu(2008, numlos, /uniform,/double)
   
   iaxis = 1+floor(rand1*3.0)
   xpos  = boxsize*rand2
   ypos  = boxsize*rand3
   zpos  = boxsize*rand4

endelse



;; Write the data
openw,1,f,width=250
for i=0, numlos-1 do begin
   printf,1,iaxis[i],xpos[i],ypos[i],zpos[i]
endfor
close,1


end



;;-------------------------------------------------------------------
