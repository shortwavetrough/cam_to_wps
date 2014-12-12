; Purpose:    Read in daily averaged CLM/POP/CIC data and copy
;             those values to create 6-hr data. Also read in 
;             invariant fields and copy to create 6-hr data.
;             Then append invariant fields to cam-6hr file, and 
;             extract individual times for CAM_netcdf_to_GRIB.ncl.
;
; Programmer: Brandon Fisel
;             Department of Geological and Atmospheric Sciences
;             Iowa State University
;
; Date:       10/27/2014
;
; Use:        ncl < copy_times.ncl
;
; Input:      CLM/POP/CIC + invariant fields (NetCDF)
;
; Output:     6-hr data (NetCDF)
;
; References: NCL reference webpage: http://www.ncl.ucar.edu/
;
; Load NCL library files
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
;
begin

; Main input parameters:
;************************************************************************
  ;directory containing data
  cdir = "$WORKDIR/intermediate/"
  ;directory to write 6-hr files
  odir = "$WORKDIR/intermediate/"
  ;directory to write individual files
  inddir = "$WORKDIR/ind/"
  ;case name to be part of the individual file filenames
  casename = "CCSM4_historical"
  ;loop over years
  do year=1979, 1989
    ;file names of each model output by ext_dates.sh, not including .nc
    cam = "cam_"+year
    clm = "clm_"+year
    pop = "pop_"+year
    cic = "cic_"+year
    inv1 = "invariant/landfrac" ;landfrac invariant file
    inv2 = "invariant/phis" ;phis invariant file
  
; Open data
;************************************************************************
    fcam = addfile(cam+".nc","r")
    fclm = addfile(clm+".nc","r")
    fpop = addfile(pop+".nc","r")
    fcic = addfile(cic+".nc","r")
    finv1 = addfile(inv1+".nc","r")
    finv2 = addfile(inv2+".nc","r")

; Get variables and process
;************************************************************************
    ;CAM time
    ftime = doubletofloat(fcam->time) ;float
    ftime!0 = "time"
    dtime = fcam->time ;double

    ;CLM
    SOILWATER_10CM_in = fclm->SOILWATER_10CM
    TSOI_10CM_in = fclm->TSOI_10CM
    lat = fclm->lat
    lon = fclm->lon

    ;Create new variables with same number of CAM times
    SOILWATER_10CM = new((/dimsizes(ftime),192,288/),typeof(SOILWATER_10CM_in))
    TSOI_10CM = new((/dimsizes(ftime),192,288/),typeof(TSOI_10CM_in))
    ;Copy attributes, and set coord names
    copy_VarAtts(SOILWATER_10CM_in,SOILWATER_10CM)
    copy_VarAtts(TSOI_10CM_in,TSOI_10CM)
    SOILWATER_10CM!0 = "time"
    SOILWATER_10CM!1 = "lat"
    SOILWATER_10CM!2 = "lon"
    TSOI_10CM!0 = "time"
    TSOI_10CM!1 = "lat"
    TSOI_10CM!2 = "lon"
    do n=0, dimsizes(ftime) - 1
      SOILWATER_10CM(n,:,:) = SOILWATER_10CM_in(n/4,:,:)
      TSOI_10CM(n,:,:) = TSOI_10CM_in(n/4,:,:)
    end do
    SOILWATER_10CM&time = ftime
    TSOI_10CM&time = ftime

    ;write CLM file and clean up
    system("/bin/rm -f "+odir+"clm_6hrly.nc")

    ;open a new file
    ncdf = addfile(odir+"clm_6hrly.nc","c")
    filedimdef(ncdf,"time",-1,True)

    ;output variables
    ncdf->SOILWATER_10CM = SOILWATER_10CM
    ncdf->TSOI_10CM = TSOI_10CM

    delete([/SOILWATER_10CM_in,TSOI_10CM_in,lat,lon,SOILWATER_10CM,\
             TSOI_10CM,ncdf,fclm/])

    ;POP
    tos_in = fpop->tos
    lat = fpop->lat
    lon = fpop->lon

    ;Create new variables with same number of CAM times
    tos = new((/dimsizes(dtime),384,320/),typeof(tos_in))
    ;Copy attributes, and set coord names
    copy_VarAtts(tos_in,tos)
    tos!0 = "time"
    tos!1 = "nlat"
    tos!2 = "nlon"
    do n=0, dimsizes(dtime) - 1
      tos(n,:,:) = tos_in(n/4,:,:)
    end do
    tos&time = dtime

    ;write POP file and clean up
    system("/bin/rm -f "+odir+"pop_6hrly.nc")

    ;open a new file
    ncdf = addfile(odir+"pop_6hrly.nc","c")
    filedimdef(ncdf,"time",-1,True)

    ;output variables
    ncdf->lat = lat
    ncdf->lon = lon
    ncdf->tos = tos

    delete([/tos_in,lat,lon,tos,ncdf,fpop/])

    ;CIC
    aice_d_in = fcic->aice_d
    lat = fcic->lat
    lon = fcic->lon

    ;Create new variables with same number of CAM times
    aice_d = new((/dimsizes(ftime),384,320/),typeof(aice_d_in))
    ;Copy attributes, and set coord names
    copy_VarAtts(aice_d_in,aice_d)
    aice_d!0 = "time"
    aice_d!1 = "nj"
    aice_d!2 = "ni"
    do n=0, dimsizes(ftime) - 1
      aice_d(n,:,:) = aice_d_in(n/4,:,:)
    end do

    aice_d&time = ftime

    ;write CIC file and clean up
    system("/bin/rm -f "+odir+"cic_6hrly.nc")

    ;open a new file
    ncdf = addfile(odir+"cic_6hrly.nc","c")
    filedimdef(ncdf,"time",-1,True)

    ;output variables
    ncdf->lat = lat
    ncdf->lon = lon
    ncdf->aice_d = aice_d

    delete([/aice_d_in,lat,lon,aice_d,ncdf,fcic/])

    ;Invariant
    landfrac_in = finv1->LANDFRAC
    phis_in = finv2->PHIS
    lat = finv1->lat
    lon = finv1->lon

    ;Create new variables with same number of CAM times
    LANDFRAC = new((/dimsizes(dtime),192,288/),typeof(landfrac_in))
    PHIS = new((/dimsizes(dtime),192,288/),typeof(phis_in))
    ;Copy attributes, and set coord names
    copy_VarAtts(landfrac_in,LANDFRAC)
    copy_VarAtts(phis_in,PHIS)
    LANDFRAC!0 = "time"
    LANDFRAC!1 = "lat"
    LANDFRAC!2 = "lon"
    PHIS!0 = "time"
    PHIS!1 = "lat"
    PHIS!2 = "lon"
    do n=0, dimsizes(dtime) - 1
      LANDFRAC(n,:,:) = landfrac_in(0,:,:)
      PHIS(n,:,:) = phis_in(0,:,:)
    end do

    LANDFRAC&time = dtime
    PHIS&time = dtime

    ;write Invariant file and clean up
    system("/bin/rm -f "+odir+"inv_6hrly.nc")

    ;open a new file
    ncdf = addfile(odir+"inv_6hrly.nc","c")
    filedimdef(ncdf,"time",-1,True)

    ;output variables
    ncdf->lat = lat
    ncdf->lon = lon
    ncdf->LANDFRAC = LANDFRAC
    ncdf->PHIS = PHIS

    delete([/landfrac_in,phis_in,lat,lon,LANDFRAC,PHIS,ncdf,finv1,\
             finv2,ftime,dtime,fcam/])

    ;Append invariant fields to CAM 6-hrly file
    system("ncks -A -v LANDFRAC,PHIS inv_6hrly.nc cam_"+year+".nc")

    ;Write times to individual files
    ;Open CAM file to get dates needed to create file strings for writing
    ;individual file names
    fcam = addfile(cam+".nc","r")
    dates = fcam->date
    delete(fcam)
    hours = (/"00","06","12","18"/)
    ;counter to determine hour
    c1 = 0

    ;Write individual files
    do n=0, dimsizes(dates) - 1
      if(c1.lt.4.) then
        HR = hours(c1)
      else
        c1 = 0
        HR = hours(0)
      end if
      int2st_date = tostring(dates(n))
      st_date = int2st_date+HR
      system("ncks -d time,"+n+" cam_"+year+".nc "+inddir+"cam_"+casename+"_"+st_date+".nc")
      system("ncks -d time,"+n+" clm_6hrly.nc "+inddir+"clm_"+casename+"_"+st_date+".nc")
      system("ncks -d time,"+n+" pop_6hrly.nc "+inddir+"pop_"+casename+"_"+st_date+".nc")
      system("ncks -d time,"+n+" cic_6hrly.nc "+inddir+"cic_"+casename+"_"+st_date+".nc")
      c1 = c1+1
    end do

    delete([/dates,hours,c1,HR,int2st_date,st_date/])
  end do

end
