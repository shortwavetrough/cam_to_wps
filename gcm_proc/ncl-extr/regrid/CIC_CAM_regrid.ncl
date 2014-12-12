;----------------------------------------------------------------------
; ESMF_template.ncl
;
; This is a template file for use with ESMF regridding.
;
; There are four main steps to regrid data from one grid to another:
;
;   1. Generate a description file (SCRIP or ESMF) for the source grid.
;
;   2. Generate a description file (SCRIP or ESMF) for the destination
;      grid.
;
;   3. Generate the weights file, using the source and destination
;      files created in #1 and #2.
;
;   4. Apply the weights to the data you want to regrid.
;
;----------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
;---Input file that contains lat/lon grid of source grid
    srcFileName = "cic_CCSM4_historical_197909-197912-1979090100.nc"

;---Input file that contains lat/lon grid of destination grid, if applicable.
    dstFileName = "cam_CCSM4_historical_197909-197912-1979090100.nc"

;
; Output (and input) files for grid descriptions and weight file
; These can be whatever names you want.
;
    srcGridName = "src_SCRIP.CIC.nc"     ; "src_ESMF.nc"
    dstGridName = "dst_SCRIP.nc"     ; "dst_ESMF.nc"
    wgtFileName = "WRF_2_Rect.CIC.nc"   ; this is just an example

;---Set to True if you want to skip any of these steps
    SKIP_SRC_GEN = False
    SKIP_DST_GEN = False
    SKIP_WGT_GEN = False

;----------------------------------------------------------------------
; 1. Generate a description file (SCRIP or ESMF) for the source grid.
;----------------------------------------------------------------------
    src_file = addfile(srcFileName,"r")
    x        = src_file->aice_d(0,:,:)            ; You may need to subscript this.
    lat_src      = src_file->lat(:,:)
    lon_src      = src_file->lon(:,:)

    x@lat2d = lat_src
    x@lon2d = lon_src

    if(.not.SKIP_SRC_GEN) then
      Opt                = True
      Opt@ForceOverWrite = True
      Opt@PrintTimings   = True
      Opt@Title          = "WRF Grid"   ; An example

;---If we have 2D lat/lon arrays.
     curvilinear_to_SCRIP(srcGridName,src_file->lat,src_file->lon,Opt)

;---Clean up
      delete(Opt)
   end if

;----------------------------------------------------------------------
; 2. Generate a description file (SCRIP or ESMF) for the destination
;    grid.
;----------------------------------------------------------------------

;---If applicable, open file that contains destination grid.
    dst_file = addfile(dstFileName,"r")

    if(.not.SKIP_DST_GEN) then
      Opt                = True
     ; Opt@LLCorner       = (/ -60.d,   0.d/)  ; default is (/-90,-180/)
     ; Opt@URCorner       = (/  60.d, 355.d/)  ; default is (/ 90, 180/)
      Opt@ForceOverWrite = True
      Opt@PrintTimings   = True


;---If we have 1D lat/lon coordinate arrays.
     rectilinear_to_SCRIP(dstGridName,dst_file->lat,dst_file->lon,Opt)

;---If we have 2D lat/lon arrays.
;     curvilinear_to_SCRIP(dstGridName,dst_file->lat,dst_file->lon,Opt)

;---Clean up
      delete(Opt)
    end if

;----------------------------------------------------------------------
; 3. Generate the weights file, using the source and destination
;    files created in #1 and #2.
;----------------------------------------------------------------------
    if(.not.SKIP_WGT_GEN) then
      Opt                      = True
      Opt@InterpMethod         = "patch"  ; "bilinear","patch", "conserve"
      Opt@ForceOverWrite       = True
      Opt@PrintTimings         = True

      ESMF_regrid_gen_weights(srcGridName,dstGridName,wgtFileName,Opt)

;---Clean up
      delete(Opt)
    end if

;----------------------------------------------------------------------
; 4. Apply the weights to the data you want to regrid.
;----------------------------------------------------------------------
    Opt                = True
;   Opt@Debug          = True
    Opt@PrintTimings   = True

    x_regrid = ESMF_regrid_with_weights(x,wgtFileName,Opt)

;
; Add attributes and coordinate arrays for plotting.
; This section will have to be modified depending
; on your destination grid.
; 
    copy_VarAtts(x,x_regrid)    ; Copy attributes

;---Output grid is rectilinear
    dstlat       = retrieve_SCRIP_lat(dstGridName)
    dstlon       = retrieve_SCRIP_lon(dstGridName)
    dstlat@units = "degrees_north"
    dstlon@units = "degrees_east"

    x_regrid!0   = "lat"
    x_regrid!1   = "lon"
    x_regrid&lat = dstlat(:,0)    ; This is a rectilinear grid, so
    x_regrid&lon = dstlon(0,:)    ; we only need a 1D sub-selection.
 
    printVarSummary(x)
    printVarSummary(x_regrid)
 
;---Output grid is curvilinear
;    dstlat       = retrieve_SCRIP_lat(dstGridName)
;    dstlon       = retrieve_SCRIP_lon(dstGridName)
;    dstlat@units = "degrees_north"
;    dstlon@units = "degrees_east"

;----------------------------------------------------------------------
; Step 4
;     Plot the original and regridded data.
;----------------------------------------------------------------------
    wks = gsn_open_wks("pdf","./cic")             ; ESMF.ps
    gsn_define_colormap(wks,"gui_default")      ; choose colormap

    res                     = True              ; Plot mods desired.
    res@gsnDraw             = False             ; We will panel later.
    res@gsnFrame            = False

    res@gsnMaximize         = True              ; Maximize plot
    res@gsnSpreadColors     = True              ; use full range of colormap
    
    res@mpMaxLatF           =   90              ; Set accordingly
    res@mpMinLatF           =  -90
    res@mpMinLonF           = -180
    res@mpMaxLonF           =  180

    res@cnFillOn            = True              ; color plot desired
    res@cnLinesOn           = False             ; turn off contour lines
    res@cnLineLabelsOn      = False             ; turn off contour lines

;    res@cnLevelSelectionMode =  "ManualLevels"   ; manual levels
;    res@cnMinLevelValF       = 4                 ; min level
;    res@cnMaxLevelValF       = 32                ; max level
;    res@cnLevelSpacingF      =  2                ; interval

    res@lbLabelBarOn        = False              ; Labelbar will be in panel

;---Plot data on original grid
    res@gsnAddCyclic = True

    dims = tostring(dimsizes(x))
    res@tiMainString = "Original data (" + str_join(dims," x ") + ")"
    plot_orig = gsn_csm_contour_map(wks,x,res)     
    delete(dims)

;---Plot data on new grid
    res@gsnAddCyclic = True

    dims = tostring(dimsizes(x_regrid))
    res@tiMainString = "Regridded data (" + str_join(dims," x ") + ")"
    plot_regrid = gsn_csm_contour_map(wks,x_regrid,res)     

;---Resources for paneling
    pres                  = True
    pres@gsnMaximize      = True
    pres@gsnPanelLabelBar = True

    gsn_panel(wks,(/plot_orig,plot_regrid/),(/2,1/),pres)
end

