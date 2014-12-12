; Purpose:    Read in CAM netcdf individual times, calculate needed
;             variables, interpolate 3D variables to pressure levels, 
;             regrid data to CAM grid, and write data to GRIB.
;             This script was adopted from Ruby Leung's Fortran 
;             program.
;
; Programmer: Brandon Fisel
;             Department of Geological and Atmospheric Sciences
;             Iowa State University
;
; Date:       10/27/2014
;
; Use:        ncl < CAM_netcdf_to_GRIB.ncl
;
; Input:      CAM (NetCDF)
;
; Output:     Ungrid-ready GRIB files
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
  ;directory containing individual time step data
  cdir = "$WORKDIR/ind/"
  ;directory to contain grib files
  odir = "$WORKDIR/output/"
  ;file names of each model, not including the date
  mod1 = "cam_CCSM4_historical_"
  mod2 = "clm_CCSM4_historical_"
  mod3 = "cic_CCSM4_historical_"
  mod4 = "pop_CCSM4_historical_"
  ;start date
  syr = st_yr
  smo = "st_mo"
  ;Change the hours below to reflect both the first filename hour, and
  ;the interval the data is available
  shr = (/"00","06","12","18"/)
  ;hour counter, do not change
  hr_c = 0.
  ;debug?
  DIAG = False
  ;Regrid weights file name
  ;Note: all regrid files are in ./regrid directory
  wgtFileName1 = "regrid/WRF_2_Rect.CIC.nc"
  wgtFileName2 = "regrid/WRF_2_Rect.POP.nc"
  ;pressure levels to interpolate from CAM sigma levels
  pnew = (/1000.0,975.0,950.0,925.0,900.0,875.0,850.0,\
            825.0,800.0,775.0,750.0,700.0,650.0,600.0,\
            550.0,500.0,450.0,400.0,350.0,300.0,250.0,\
            225.0,200.0,175.0,150.0,125.0,100.0,70.0,\
             50.0,30.0,20.0,10.0,7.0,5.0,3.0,2.0,1.0/)
  
; Open data
;************************************************************************
  files = systemfunc("ls "+cdir+mod1+syr+smo+"*.nc")
  netcdf_cam_filenames = addfiles(files,"r")
  ;get the number of files for do loop below
  numfil = dimsizes(files) - 1
  delete(files)
  files = systemfunc("ls "+cdir+mod2+syr+smo+"*.nc")
  netcdf_clm_filenames = addfiles(files,"r")
  delete(files)
  files = systemfunc("ls "+cdir+mod3+syr+smo+"*.nc")
  netcdf_cic_filenames = addfiles(files,"r")
  delete(files)
  files = systemfunc("ls "+cdir+mod4+syr+smo+"*.nc")
  netcdf_pop_filenames = addfiles(files,"r")
  delete(files)

; Begin loop over daily files in the month
;************************************************************************
  do f=0, numfil
  netcdf_cam_filename = netcdf_cam_filenames[f]
  netcdf_clm_filename = netcdf_clm_filenames[f]
  netcdf_cic_filename = netcdf_cic_filenames[f]
  netcdf_pop_filename = netcdf_pop_filenames[f]
  ;counter to determine hour
  if (hr_c.eq.0.) then
    HR = shr(0)
  else if (hr_c.eq.1.) then
    HR = shr(1)
  else if (hr_c.eq.2.) then
    HR = shr(2)
  else if (hr_c.eq.3.) then
    HR = shr(3)
  else if (hr_c.eq.4.) then
    hr_c = 0.
    HR = shr(0)
  else
    print("Hour counter is off. Exiting...")
    status_exit(42)
  end if
  end if
  end if
  end if
  end if
; Get variables
;************************************************************************
  ;CAM
  lat = netcdf_cam_filename->lat
  lon = netcdf_cam_filename->lon
  lev = netcdf_cam_filename->lev
  time = netcdf_cam_filename->time
  date = netcdf_cam_filename->date
  date_to_str = tostring(date)
  dtime_map = (/4,2,2/)
  st_dtime = str_split_by_length(date_to_str, dtime_map)
  YYYY = st_dtime(0)
  MM = st_dtime(1)
  DD = st_dtime(2)
  nlon = netcdf_cam_filename->nlon
  T = netcdf_cam_filename->T(0,:,:,:)
  Q = netcdf_cam_filename->Q(0,:,:,:)
  U = netcdf_cam_filename->U(0,:,:,:)
  V = netcdf_cam_filename->V(0,:,:,:)
  Z3 = netcdf_cam_filename->Z3(0,:,:,:)
  PS = netcdf_cam_filename->PS(0,:,:)
  PSL = netcdf_cam_filename->PSL(0,:,:)
  LANDFRAC = netcdf_cam_filename->LANDFRAC(0,:,:)
  PHIS = netcdf_cam_filename->PHIS(0,:,:)
  P0 = netcdf_cam_filename->P0
  hyam = netcdf_cam_filename->hyam
  hybm = netcdf_cam_filename->hybm

  ;CLM
  SOILWATER_10CM = netcdf_clm_filename->SOILWATER_10CM(0,:,:)
  TSOI_10CM = netcdf_clm_filename->TSOI_10CM(0,:,:)

  ;CIC
  aice_d = netcdf_cic_filename->aice_d(0,:,:)

  ;POP
  TOS = netcdf_pop_filename->tos(0,:,:)

  delete([/netcdf_cam_filename,netcdf_clm_filename,\
           netcdf_cic_filename,netcdf_pop_filename/])
; Create new arrays to hold calculated variables
; and set units and long name
;************************************************************************
  lat@units = "degrees_north"
  lat@long_name = "latitude"
  lon@units = "degrees_east"
  lon@long_name = "longitude"

  i3d_TT = new((/37,192,288/),"float")
  i3d_TT!0 = "lev_p"
  i3d_TT!1 = "lat"
  i3d_TT!2 = "lon"
  i3d_TT@units = "K"
  i3d_TT@long_name = "CAM temperature"
  i3d_RH = new((/37,192,288/),"float")
  i3d_RH@units = "%"
  i3d_RH@long_name = "RH computed from CAM specific humidity"
  i3d_UU = new((/37,192,288/),"float")
  i3d_UU!0 = "lev_p"
  i3d_UU!1 = "lat"
  i3d_UU!2 = "lon"
  i3d_UU@units = "m/s"
  i3d_UU@long_name = "CAM U"
  i3d_VV = new((/37,192,288/),"float")
  i3d_VV!0 = "lev_p"
  i3d_VV!1 = "lat"
  i3d_VV!2 = "lon"
  i3d_VV@units = "m/s"
  i3d_VV@long_name = "CAM V"
  i3d_GHT = new((/37,192,288/),"float")
  i3d_GHT!0 = "lev_p"
  i3d_GHT!1 = "lat"
  i3d_GHT!2 = "lon"
  i3d_GHT@units = "m"
  i3d_GHT@long_name = "Geopotential Height (above sea level)"

  i2d_PSFC = new((/192,288/),"float")
  i2d_PSFC!0 = "lat"
  i2d_PSFC!1 = "lon"
  i2d_PSFC@units = "Pa"
  i2d_PSFC@long_name = "Surface Pressure"
  i2d_PMSL = new((/192,288/),"float")
  i2d_PMSL!0 = "lat"
  i2d_PMSL!1 = "lon"
  i2d_PMSL@units = "Pa"
  i2d_PMSL@long_name = "Sea level pressure"
  i2d_landsea = new((/192,288/),"float")
  i2d_landsea!0 = "lat"
  i2d_landsea!1 = "lon"
  i2d_landsea@units = "percent"
  i2d_landsea@long_name = "land mask from CLM"
  i2d_SKINTEMP = new((/192,288/),"float")
  i2d_SKINTEMP!0 = "lat"
  i2d_SKINTEMP!1 = "lon"
  i2d_SKINTEMP@units = "K"
  i2d_SKINTEMP@long_name = "regridded SST from POP to CAM"
  i2d_TT = new((/192,288/),"float")
  i2d_TT@units = "K"
  i2d_TT@long_name = "VerticalExtrapolate to 2m from CAM T at lev=26"
  i2d_RH = new((/192,288/),"float")
  i2d_RH@units = "%"
  i2d_RH@long_name = "compute RH from CAM SPECHUMD at lev=26"
  i2d_UU = new((/192,288/),"float")
  i2d_UU!0 = "lat"
  i2d_UU!1 = "lon"
  i2d_UU@units = "m/s"
  i2d_UU@long_name = "VerticalExtrapolate to 10m from CAM U at lev26"
  i2d_VV = new((/192,288/),"float")
  i2d_VV!0 = "lat"
  i2d_VV!1 = "lon"
  i2d_VV@units = "m/s"
  i2d_VV@long_name = "VerticalExtrapolate to 10m from CAM V at lev26"
  i2d_SEAICE = new((/192,288/),"float")
  i2d_SEAICE!0 = "lat"
  i2d_SEAICE!1 = "lon"
  i2d_SEAICE@units = "percent"
  i2d_SEAICE@long_name = "CIC ice area (aggregate)"

  SM000010 = new((/192,288/),"float")
  SM000010!0 = "lat"
  SM000010!1 = "lon"
  SM000010@units = "mm3/mm3"
  SM000010@long_name = "volumetric soil water"
  SM010040 = new((/192,288/),"float")
  SM010040!0 = "lat"
  SM010040!1 = "lon"
  SM010040@units = "mm3/mm3"
  SM010040@long_name = "volumetric soil water"
  SM040100 = new((/192,288/),"float")
  SM040100!0 = "lat"
  SM040100!1 = "lon"
  SM040100@units = "mm3/mm3"
  SM040100@long_name = "volumetric soil water"
  SM100200 = new((/192,288/),"float")
  SM100200!0 = "lat"
  SM100200!1 = "lon"
  SM100200@units = "mm3/mm3"
  SM100200@long_name = "volumetric soil water"
  ST000010 = new((/192,288/),"float")
  ST000010!0 = "lat"
  ST000010!1 = "lon"
  ST000010@units = "K"
  ST000010@long_name = "soil temperature"
  ST010040 = new((/192,288/),"float")
  ST010040!0 = "lat"
  ST010040!1 = "lon"
  ST010040@units = "K"
  ST010040@long_name = "soil temperature"
  ST040100 = new((/192,288/),"float")
  ST040100!0 = "lat"
  ST040100!1 = "lon"
  ST040100@units = "K"
  ST040100@long_name = "soil temperature"
  ST100200 = new((/192,288/),"float")
  ST100200!0 = "lat"
  ST100200!1 = "lon"
  ST100200@units = "K"
  ST100200@long_name = "soil temperature"

; Regrid CIC/POP to CAM grid
;************************************************************************
; This uses the new ECMF regridding functions. There are four main steps
; to regrid data from one grid to another:
; 1) Generate a description file (SCRIP or ESMF) for the source grid,
; 2) Generate a description file (SCRIP or ESMF) for the destination 
;    grid,
; 3) Generate the weights file, using the source and destination files 
;    created in (1) and (2), and
; 4) Apply the weights to the data you want to regrid.
; Note: These steps are completed using NCL files in /regrid directory. 
; Below only uses the output from those scripts in the /regrid directory.
;************************************************************************
  Opt = True
  Opt@Debug = True
  Opt@PrintTimings = True

  cic_regrid = ESMF_regrid_with_weights(aice_d,wgtFileName1,Opt)
  pop_regrid = ESMF_regrid_with_weights(TOS,wgtFileName2,Opt)

  ;Add attributes and coordinate arrays for plotting and will have to be
  ;modified depending on your destination grid.
  copy_VarAtts(aice_d,cic_regrid)
  copy_VarAtts(TOS,pop_regrid)

  ;Output grid is rectilinear
  dstlat = retrieve_SCRIP_lat("regrid/dst_SCRIP.nc")
  dstlon = retrieve_SCRIP_lon("regrid/dst_SCRIP.nc")
  dstlat@units = "degrees_north"
  dstlon@units = "degrees_east"

  cic_regrid!0 = "lat"
  cic_regrid!1 = "lon"
  ;This is a rectilinear grid, so only need a 1D sub-selection
  cic_regrid&lat = dstlat(:,0)
  cic_regrid&lon = dstlon(0,:)

  pop_regrid!0 = "lat"
  pop_regrid!1 = "lon"
  ;This is a rectilinear grid, so only need a 1D sub-selection
  pop_regrid&lat = dstlat(:,0)
  pop_regrid&lon = dstlon(0,:)

  delete([/aice_d, TOS, Opt, dstlat, dstlon/])
; Calculate vaiables needed by WPS, and fix other variables
;************************************************************************
  i2d_PSFC = PS
  i2d_PMSL = PSL

  print("Fixing landfrac...")
  i2d_landsea = where(LANDFRAC.lt.0.5,0,1)
  if (DIAG) then
    print("land mask after fixing: "+i2d_landsea(:,64))
  end if
  print("Done fixing landfrac to be only 0 or 1")

  print("Fixing seaice to be fractional...")
  i2d_SEAICE = cic_regrid*0.01
  if (DIAG) then
    print("seaice after fixing: "+i2d_SEAICE(:,64))
  end if
  print("Done fixing seaice to be fractional")

  print("Fixing SKINTEMP to be in Kelvin...")
  i2d_SKINTEMP = pop_regrid+273.15
  if (DIAG) then
    print("skintemp after fixing: "+i2d_SKINTEMP(:,64))
  end if
  print("Done fixing skintemp to be in Kelvin")

  print("Calculating topo...")
  topo = PHIS/9.81
  if (DIAG) then
    print("topo calculated: "+topo(:,64))
  end if
  print("Done calculating topo")

  print("Calculating 2-d surface fields...")
  tbot = T(25,:,:)
  ubot = U(25,:,:)
  vbot = V(25,:,:)
  z3bot = Z3(25,:,:)

  i2d_UU = 10.*ubot/(z3bot-topo) ;U10m
  i2d_VV = 10.*vbot/(z3bot-topo) ;V10m
  t2m = 0.0065*(z3bot-(topo+2))+tbot ;T2m
  tsfc = 0.0065*(z3bot-topo)+tbot ;surface temperature, use it on land as SKINTEMP
  i2d_TT = t2m ;T2m
  i2d_TT!0 = "lat"
  i2d_TT!1 = "lon"

  RH = new(dimsizes(Q),"float")
  do k=0,25
    esat = 611.2*exp(17.67*(T(k,:,:)-273.15)/(T(k,:,:)-273.15+243.5))
    pressure = hyam(k)*P0+hybm(k)*PS
    qsat = doubletofloat(0.622*esat/(pressure-0.378*esat))
    RH(k,:,:) = (Q(k,:,:)/qsat)*100. ;RH
  end do
  i2d_RH = RH(25,:,:) ;RH-2m
  i2d_RH!0 = "lat"
  i2d_RH!1 = "lon"

  ;substitute SKINTEMP on land with tsfc, using landsea
  i2d_SKINTEMP = where(i2d_landsea.eq.1,tsfc,i2d_SKINTEMP)

  print("Fixing SKINTEMP missing values...")
  ;linearly interpolate to fill in missing values
  i2d_SKINTEMP_tmp = linmsg(i2d_SKINTEMP,0)
  delete(i2d_SKINTEMP)
  i2d_SKINTEMP = i2d_SKINTEMP_tmp
  i2d_SKINTEMP!0 = "lat"
  i2d_SKINTEMP!1 = "lon"
  delete(i2d_SKINTEMP_tmp)

  ;print some diagnostics about 2-d fields
  if (DIAG) then
    print(" --- Surface fields --- ")
    print("UU: "+i2d_UU(:,64)+", VV: "+i2d_VV(:,64)+", TT: "+i2d_TT(:,64)+\
          ", tsfc: "+tsfc(:,64)+", SKINTEMP: "+i2d_SKINTEMP(:,64))
  end if
  print("Done calculating surface fields")

  print("Setting soil variables from CLM...")
  ;set WRF soil moisture values at all levels to be equal to the first
  ;CAM soil moisture level
  ;All soil moisture should have the unit as "fraction"
  ;CCSM soil moisture is in kg water / m2
  SM000010 = SOILWATER_10CM/100.
  SM010040 = SOILWATER_10CM/100.
  SM040100 = SOILWATER_10CM/100.
  SM100200 = SOILWATER_10CM/100.

  ;set WRF soil temperature at all levels to be equal to the first CAM
  ;soil level
  ST000010 = TSOI_10CM
  ST010040 = TSOI_10CM
  ST040100 = TSOI_10CM
  ST100200 = TSOI_10CM

  ;deal with missing soil values (1.0e36)
  ;for moisture, set missing values to 0
  SM000010 = where(SM000010.gt.1.0e30,0.,SM000010)
  SM010040 = where(SM010040.gt.1.0e30,0.,SM010040)
  SM040100 = where(SM040100.gt.1.0e30,0.,SM040100)
  SM100200 = where(SM100200.gt.1.0e30,0.,SM100200)
  SM000010 = where(SM000010.gt.1.,1.,SM000010)
  SM010040 = where(SM010040.gt.1.,1.,SM010040)
  SM040100 = where(SM040100.gt.1.,1.,SM040100)
  SM100200 = where(SM100200.gt.1.,1.,SM100200)

  ;set missing values to surface temperature
  ST000010 = where(ST000010.gt.1.0e30,i2d_SKINTEMP,ST000010)
  ST010040 = where(ST010040.gt.1.0e30,i2d_SKINTEMP,ST010040)
  ST040100 = where(ST040100.gt.1.0e30,i2d_SKINTEMP,ST040100)
  ST100200 = where(ST100200.gt.1.0e30,i2d_SKINTEMP,ST100200)
  if (DIAG) then
    print(" --- Soil fields --- ")
    print("SM000010(21,121): "+SM000010(21,121)+\
          ", SM010040(21,121): "+SM010040(21,121)+\
          ", SM040100(21,121): "+SM040100(21,121)+\
          ", SM100200(21,121): "+SM100200(21,121))
    print("ST000010(21,121): "+ST000010(21,121)+\
          ", ST010040(21,121): "+ST010040(21,121)+\
          ", ST040100(21,121): "+ST040100(21,121)+\
          ", ST100200(21,121): "+ST100200(21,121))
  end if
  print("Done setting soil variables")

  print("Interpolating 3-d variables to pressure levels...")
  P0 = P0/100. ;vinth2p wants P0 in hpa
  i3d_TT = vinth2p(T,hyam,hybm,pnew,PS,1,P0,2,True)
  i3d_RH = vinth2p(RH,hyam,hybm,pnew,PS,1,P0,2,True)
  i3d_UU = vinth2p(U,hyam,hybm,pnew,PS,1,P0,2,True)
  i3d_VV = vinth2p(V,hyam,hybm,pnew,PS,1,P0,2,True)
  i3d_GHT = vinth2p(Z3,hyam,hybm,pnew,PS,1,P0,2,True)

  i3d_RH!0 = "lev_p"
  i3d_RH!1 = "lat"
  i3d_RH!2 = "lon"

  ;check RH is not negative
  c1=0.
  c2=0.
  do k=0,36
    do i=0,191
      do j=0,287
        if (i3d_RH(k,i,j).lt.0.) then
          if (c2.lt.100.) then
            print("*** RH<0: RH="+i3d_RH(k,i,j)+",k="+k+",i="+i+",j="+j)
          end if
          if (c2.eq.100.) then
            print("*** The number of negative RH points has exceeded 100.")
            print("Further output will be suppressed ***")
          end if
          i3d_RH(k,i,j) = 0.
          c2=c2+1.
        end if
        c1=c1+1.
      end do
    end do
  end do
  if (DIAG) then
    print(" --- 3d fields --- ")
    print("i3d_TT min:"+min(i3d_TT)+"i3d_TT max:"+max(i3d_TT)+\
          "i3d_TT(:,:,64):"+i3d_TT(:,:,64))
    print("i3d_RH min:"+min(i3d_RH)+"i3d_RH max:"+max(i3d_RH)+\
          "i3d_RH(:,:,64):"+i3d_RH(:,:,64))
    print("i3d_UU min:"+min(i3d_UU)+"i3d_UU max:"+max(i3d_UU)+\
          "i3d_UU(:,:,64):"+i3d_UU(:,:,64))
    print("i3d_VV min:"+min(i3d_VV)+"i3d_VV max:"+max(i3d_VV)+\
          "i3d_VV(:,:,64):"+i3d_VV(:,:,64))
    print("i3d_GHT min:"+min(i3d_GHT)+"i3d_GHT max:"+max(i3d_GHT)+\
          "i3d_GHT(:,:,64):"+i3d_GHT(:,:,64))
    print("c1(all): "+c1+", c2(neg): "+c2+" ===> "+(c2/c1)*100.+"% are negative points")
  end if
  print("Done interpolating 3-d variables to pressure levels")

  delete([/T,Q,U,V,Z3,PS,PSL,LANDFRAC,PHIS,P0,hyam,hybm,\
           SOILWATER_10CM,TSOI_10CM,cic_regrid,pop_regrid,\
           topo,tbot,ubot,vbot,z3bot,t2m,tsfc,RH,esat,pressure,\
           qsat,c1,c2/])

  print("Writing fields to NetCDF...")
; Write variables to NetCDF, then using system command to convert 
; NetCDF file to GRIB
;************************************************************************
  ;these vars need dimensions renamed here
  i3d_RH!0 = "lev_p"
  i3d_RH!1 = "lat"
  i3d_RH!2 = "lon"
  i2d_TT!0 = "lat"
  i2d_TT!1 = "lon"
  i2d_RH!0 = "lat"
  i2d_RH!1 = "lon"

  ;create directories
  system("mkdir "+odir+"netcdf/")
  system("mkdir "+odir+"grib/")
  ;remove old files
  system("/bin/rm -f "+odir+"netcdf/CAM_"+date+HR+".nc")
  system("/bin/rm -f "+odir+"grib/split_*.grib")
  system("/bin/rm -f "+odir+"grib/CAM_"+date+HR+".grib")

  ;open a new file
  ncdf = addfile(odir+"netcdf/CAM_"+date+HR+".nc","c")

  ;output variables
  ncdf->TT = i3d_TT
  ncdf->RH = i3d_RH
  ncdf->UU = i3d_UU
  ncdf->VV = i3d_VV
  ncdf->GHT = i3d_GHT
  ncdf->PSFC = i2d_PSFC
  ncdf->PMSL = i2d_PMSL
  ncdf->landsea = i2d_landsea
  ncdf->SKINTEMP = i2d_SKINTEMP
  ncdf->T2M = i2d_TT
  ncdf->RH2M = i2d_RH
  ncdf->UU10M = i2d_UU
  ncdf->VV10M = i2d_VV
  ncdf->SEAICE = i2d_SEAICE
  ncdf->SM000010 = SM000010
  ncdf->SM010040 = SM010040
  ncdf->SM040100 = SM040100
  ncdf->SM100200 = SM100200
  ncdf->ST000010 = ST000010
  ncdf->ST010040 = ST010040
  ncdf->ST040100 = ST040100
  ncdf->ST100200 = ST100200

  print("Finished writing fields to NetCDF")

  print("Creating grib file...")
; There are 4 steps to generate CAM grib files from netcdf:
; 1)Split variable names into individual files:
;   >cdo splitname <file.nc> split_
; 2)Convert the individual netcdf files to grib, and map the variables:
;   Options:
;     -setcode - GRIB1 Parameter Number
;     -setzaxis - For soil levels only; From Level1, Frome Level2
;     -setltype - Level Type
;     -setlevel - For surface fields only; From Level1
;                 e.g., 2m-T, 10m-U/V. This only makes changes to Level1 
;                 and not to Level2. Use -setzaxis to make make changes 
;                 to both levels.
;     -settabnum - Use GRIB Parameter Table #130
;   >cdo -f grb -setcode,<#> -setltype,<#> -setlevel,<#> -settabnum,130 
;      split_<var>.nc split_<var>.grib
; 3)Change GRIB times:
;   >cdo -settaxis,<YYYY-MM-DD>,<HH> split_<var>.grib split_<var>_new.grib
; 4)Merge the new GRIB files into one file:
;   >cdo merge split_*_new.grib CAM_<YYYYMMDDHH>.grib
;************************************************************************
  ;Step 1
  system("cdo splitname "+odir+"netcdf/CAM_"+date+HR+".nc "+odir+"grib/split_")

  ;Step 2
  system("cdo -f grb -setcode,1 -setltype,1 -setlevel,0 -settabnum,130 "+odir+"grib/split_PSFC.nc "+odir+"grib/split_PSFC.grib")
  system("cdo -f grb -setcode,2 -setltype,102 -setlevel,0 -settabnum,130 "+odir+"grib/split_PMSL.nc "+odir+"grib/split_PMSL.grib")
  system("cdo -f grb -setcode,81 -setltype,1 -setlevel,0 -settabnum,130 "+odir+"grib/split_landsea.nc "+odir+"grib/split_landsea.grib")
  system("cdo -f grb -setcode,11 -setltype,1 -setlevel,0 -settabnum,130 "+odir+"grib/split_SKINTEMP.nc "+odir+"grib/split_SKINTEMP.grib")
  system("cdo -f grb -setcode,11 -setltype,105 -setlevel,2 -settabnum,130 "+odir+"grib/split_T2M.nc "+odir+"grib/split_T2M.grib")
  system("cdo -f grb -setcode,52 -setltype,105 -setlevel,2 -settabnum,130 "+odir+"grib/split_RH2M.nc "+odir+"grib/split_RH2M.grib")
  system("cdo -f grb -setcode,33 -setltype,105 -setlevel,10 -settabnum,130 "+odir+"grib/split_UU10M.nc "+odir+"grib/split_UU10M.grib")
  system("cdo -f grb -setcode,34 -setltype,105 -setlevel,10 -settabnum,130 "+odir+"grib/split_VV10M.nc "+odir+"grib/split_VV10M.grib")
  system("cdo -f grb -setcode,91 -setltype,1 -setlevel,0 -settabnum,130 "+odir+"grib/split_SEAICE.nc "+odir+"grib/split_SEAICE.grib")
  system("cdo -f grb -setcode,11 -setltype,100 -settabnum,130 "+odir+"grib/split_TT.nc "+odir+"grib/split_TT.grib")
  system("cdo -f grb -setcode,52 -setltype,100 -settabnum,130 "+odir+"grib/split_RH.nc "+odir+"grib/split_RH.grib")
  system("cdo -f grb -setcode,33 -setltype,100 -settabnum,130 "+odir+"grib/split_UU.nc "+odir+"grib/split_UU.grib")
  system("cdo -f grb -setcode,34 -setltype,100 -settabnum,130 "+odir+"grib/split_VV.nc "+odir+"grib/split_VV.grib")
  system("cdo -f grb -setcode,7 -setltype,100 -settabnum,130 "+odir+"grib/split_GHT.nc "+odir+"grib/split_GHT.grib")
  system("cdo -f grb -setcode,144 -setzaxis,myzaxis0_10.dat -setltype,112 -settabnum,130 "+odir+"grib/split_SM000010.nc "+odir+"grib/split_SM000010.grib")
  system("cdo -f grb -setcode,144 -setzaxis,myzaxis10_40.dat -setltype,112 -settabnum,130 "+odir+"grib/split_SM010040.nc "+odir+"grib/split_SM010040.grib")
  system("cdo -f grb -setcode,144 -setzaxis,myzaxis40_100.dat -setltype,112 -settabnum,130 "+odir+"grib/split_SM040100.nc "+odir+"grib/split_SM040100.grib")
  system("cdo -f grb -setcode,144 -setzaxis,myzaxis100_200.dat -setltype,112 -settabnum,130 "+odir+"grib/split_SM100200.nc "+odir+"grib/split_SM100200.grib")
  system("cdo -f grb -setcode,11 -setzaxis,myzaxis0_10.dat -setltype,112 -settabnum,130 "+odir+"grib/split_ST000010.nc "+odir+"grib/split_ST000010.grib")
  system("cdo -f grb -setcode,11 -setzaxis,myzaxis10_40.dat -setltype,112 -settabnum,130 "+odir+"grib/split_ST010040.nc "+odir+"grib/split_ST010040.grib")
  system("cdo -f grb -setcode,11 -setzaxis,myzaxis40_100.dat -setltype,112 -settabnum,130 "+odir+"grib/split_ST040100.nc "+odir+"grib/split_ST040100.grib")
  system("cdo -f grb -setcode,11 -setzaxis,myzaxis100_200.dat -setltype,112 -settabnum,130 "+odir+"grib/split_ST100200.nc "+odir+"grib/split_ST100200.grib")

  ;Step 3
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_PSFC.grib "+odir+"grib/split_PSFC_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_PMSL.grib "+odir+"grib/split_PMSL_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_landsea.grib "+odir+"grib/split_landsea_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_SKINTEMP.grib "+odir+"grib/split_SKINTEMP_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_T2M.grib "+odir+"grib/split_T2M_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_RH2M.grib "+odir+"grib/split_RH2M_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_UU10M.grib "+odir+"grib/split_UU10M_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_VV10M.grib "+odir+"grib/split_VV10M_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_SEAICE.grib "+odir+"grib/split_SEAICE_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_TT.grib "+odir+"grib/split_TT_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_RH.grib "+odir+"grib/split_RH_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_UU.grib "+odir+"grib/split_UU_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_VV.grib "+odir+"grib/split_VV_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_GHT.grib "+odir+"grib/split_GHT_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_SM000010.grib "+odir+"grib/split_SM000010_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_SM010040.grib "+odir+"grib/split_SM010040_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_SM040100.grib "+odir+"grib/split_SM040100_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_SM100200.grib "+odir+"grib/split_SM100200_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_ST000010.grib "+odir+"grib/split_ST000010_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_ST010040.grib "+odir+"grib/split_ST010040_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_ST040100.grib "+odir+"grib/split_ST040100_new.grib")
  system("cdo settaxis,"+YYYY+"-"+MM+"-"+DD+","+HR+" "+odir+"grib/split_ST100200.grib "+odir+"grib/split_ST100200_new.grib")

  ;Step 4
  system("cdo merge "+odir+"grib/split_*_new.grib "+odir+"grib/CAM_"+YYYY+"-"+MM+"-"+DD+"-"+HR+".grib")

  print("Finished creating grib")
  system("/bin/rm -f "+odir+"grib/split_*")

  print("Advancing counter and cleaning up...")
  ;Advance the hour counter
  hr_c = hr_c+1.

  ;Clean up
  delete([/ncdf,i3d_TT,i3d_RH,i3d_UU,i3d_VV,i3d_GHT,i2d_PSFC,i2d_PMSL,\
           i2d_landsea,i2d_SKINTEMP,i2d_TT,i2d_RH,i2d_UU,i2d_VV,\
           i2d_SEAICE,SM000010,SM010040,SM040100,SM100200,ST000010,\
           ST010040,ST040100,ST100200/])

  ;daily loop
  end do

end
