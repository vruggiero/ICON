; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;************************************************************************
; water_budget.ncl
;
; PURPOSE:
; =========
; - calculate the water budget from the ICON triangular grid output
;   For this, the shell script water_budget.s is called (make sure that 
;   this script is located in the same directory as your NCL script). 
; - visualize ASCII output generated by the script water_budget.s
;
; AUTHOR: 
; Martin Koehler, DWD
; Daniel Reinert, DWD
;
; LAST CHANGES: 
; 2014-06-27
;****************************************************************************
begin

 ; name of file produced by the shell script water_budget.s
 filename="budget.txt"

 ; file with model level output on native grid:
 ICON_filename_ml = "????????????????????????????????????"


 ;
 ; use CDO to compute water budget
 ;
 system("./water_budget.s ./????????????????????????????? >budget.txt")

 ; Open model level output file
 ; (only required for reading the time axis)
 File_ml = addfile( ICON_filename_ml, "r" )

 ; read time axis from NetCDF file
 if (isfilevar(File_ml, "time")) then
   time = File_ml->time(:)   ; dims: (time)
   print("reading time complete" )
   ntimes = dimsizes(time)                   ; number of time slices
 else
   print("WARNING: time not found in input file")
 end if


 ; read ASCII data from output of water budget script
 ;
 data = asciiread(filename,-1, "string")



 ; read dingle columns to separate arrays
 ; (date and time columns are skipped) 
 offset = 4               ; number of lines to skip in ASCII file
 ;
 lhf      = stringtofloat(str_get_field(data(offset:offset+ntimes-1),  3," "))
 tot_prec = stringtofloat(str_get_field(data(offset:offset+ntimes-1),  4," "))
 tqv      = stringtofloat(str_get_field(data(offset:offset+ntimes-1),  5," "))
 tqc      = stringtofloat(str_get_field(data(offset:offset+ntimes-1),  6," "))
 tqi      = stringtofloat(str_get_field(data(offset:offset+ntimes-1),  7," "))
 tqr      = stringtofloat(str_get_field(data(offset:offset+ntimes-1),  8," "))
 tqs      = stringtofloat(str_get_field(data(offset:offset+ntimes-1),  9," "))
 pme      = stringtofloat(str_get_field(data(offset:offset+ntimes-1), 10," "))
 dqt      = stringtofloat(str_get_field(data(offset:offset+ntimes-1), 11," "))
 budget   = stringtofloat(str_get_field(data(offset:offset+ntimes-1), 12," "))

 ; add name to cordinate 0
 lhf!0       = "time"
 tot_prec!0  = "time"
 tqv!0       = "time"
 tqc!0       = "time"
 tqi!0       = "time"
 tqr!0       = "time"
 tqs!0       = "time"
 pme!0       = "time"
 dqt!0       = "time"
 budget!0    = "time"


 ; assign time array to coordinate axis 0
 lhf&time      = time
 tot_prec&time = time
 tqv&time      = time
 tqc&time      = time
 tqi&time      = time
 tqr&time      = time
 tqs&time      = time
 pme&time      = time
 dqt&time      = time
 budget&time   = time



;*********************************************************************
;  plotting
;*********************************************************************

 wks   = gsn_open_wks ("ps","water_budget")      ; open workstation

 res                  = True                     ; plot mods desired
 res@gsnDraw          = False           
 res@gsnFrame         = False   
 res@tiMainString     = "Basic XY plot"          ; add title

 res@vpHeightF= 0.47                             ; change aspect ratio of plot
 res@vpWidthF = 0.75             
 res@xyLineThicknessF = 1.5                      ; change line thickness

 res@tiMainString     = ""
 res@tiYAxisString    = "column integrated water vapour [mm]"    ; y-axis label    
 res@tiXAxisString    = "Forecast time [min]"    ; y-axis label   


;**********************************************
; time series of column integrated water vapour
;**********************************************
 plot_tqv = gsn_csm_xy (wks,tqv&time,tqv,res) ; create plot

 draw(plot_tqv)
 frame(wks)



;**********************************************
; time series of column integrated condensates
;**********************************************
 data_tq = new((/4,dimsizes(tqc&time)/),float)


 res@xyLineThicknesses = (/1.6,1.6,1.6,1.6/)
 res@xyLineColors      = (/"blue","red","darkgreen","black"/)
 res@xyDashPatterns    = (/0, 1, 7, 12/)
 res@trYMaxF           =  0.065                 ; max value on y-axis

 res@tiYAxisString    = "column integrated condensate [mm]"    ; y-axis label      
 res@tiXAxisString    = "Forecast time [min]"    ; y-axis label   

 data_tq(0,:) = tqc
 data_tq(1,:) = tqi
 data_tq(2,:) = tqr
 data_tq(3,:) = tqs
 plot_tq = gsn_csm_xy (wks,data_tq&time,data_tq,res) ; create plot

;
; legend resources
;
 lgres                    = True
 lgres@lgLineColors       = (/"blue","red","darkgreen","black"/)    ; colors for legend lines
 lgres@vpWidthF           = 0.2                   ; width of legend (NDC)
 lgres@vpHeightF          = 0.13                   ; height of legend (NDC)
 lgres@lgPerimOn          = False                 ; turn off perimeter
 lgres@lgDashIndexes      = (/0,1,7,12/)              ; dash indexes
 lgres@lgLineThicknessF   = 1.9
 lgres@lgLineLabelStrings = (/"","","",""/)             ; no labels in lines

 gsn_legend_ndc(wks,4,(/"tqc","tqi","tqr","tqs"/),0.75,0.80,lgres)


 draw(plot_tq)
 frame(wks)


;**********************************************
; P-E and dqt
;**********************************************

 data_bud = new((/2,dimsizes(pme&time)/),float)

 ; compute hourly delta
 do n=1, dimsizes(pme&time)-1
  pme(n) = pme(n)/tofloat((time(n)-time(n-1))/60.)
  dqt(n) = dqt(n)/tofloat((time(n)-time(n-1))/60.)
 end do

 data_bud(0,:) = pme
 data_bud(1,:) = dqt

 res@trYMaxF          =  0.025                 ; max value on y-axis
 res@trYMinF          = -0.025                 ; min value on y-axis
 res@tiYAxisString    = "Deviation [mm/h]"     ; y-axis label      
 res@tiXAxisString    = "Forecast time [min]"  ; y-axis label   
 res@gsnYRefLine      = 0.0                    ; create a reference line   
 res@gsnCenterString  = "(Precip - Evap) and total moisture change (dqt) per hour"

 plot_bud  = gsn_csm_xy (wks,data_bud&time,data_bud,res) ; create plot

; legend resources
;
 delete(lgres)
 lgres                    = True
 lgres@lgLineColors       = (/"blue","red"/)    ; colors for legend lines
 lgres@vpWidthF           = 0.17                   ; width of legend (NDC)
 lgres@vpHeightF          = 0.10                   ; height of legend (NDC)
 lgres@lgPerimOn          = False                 ; turn off perimeter
 lgres@lgDashIndexes      = (/0,1/)              ; dash indexes
 lgres@lgLineThicknessF   = 1.9
 lgres@lgLineLabelStrings = (/"",""/)             ; no labels in lines

 gsn_legend_ndc(wks,2,(/"P-E","dqt"/),0.72,0.80,lgres)

 draw(plot_bud)
 frame(wks)

;**********************************************
; budget (residuum)
;**********************************************

data_res =  new((/1,dimsizes(budget&time)/),float)

data_res = budget        ; just to get metadata

do n=0, dimsizes(budget&time)-1
  data_res(0,n) =  sum(budget(0:n))
end do

 res@trYMaxF          =  0.03                 ; max value on y-axis
 res@trYMinF          = -0.03                 ; min value on y-axis
 res@tiYAxisString    = "Residuum [mm]"    ; y-axis label      
 res@tiXAxisString    = "Forecast time [min]"    ; y-axis label   
 res@gsnYRefLine      = 0.0             ; create a reference line   
 res@gsnCenterString  = "Integrated water budget residuum"

 plot_res = gsn_csm_xy (wks,data_res&time,data_res,res) ; create plot

 draw(plot_res)
 frame(wks)
end
