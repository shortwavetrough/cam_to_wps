#!/bin/bash

#Pre-process GCM data to be later ingested into copy_times.ncl by
#creating yearly files containing only necessary variables.

#Directory to store yearly files
cd $WORKDIR/intermediate/
#Directory to store temporary files
mkdir tmp/

for year in {1979..1989} #loop through all years
do
  yearb=$(($year-1))

  #CLM data
  ncks -v SOILWATER_10CM,TSOI_10CM,topo -d time,364 "../data/b40.20th.track1.1deg.012.clm2.h1."$yearb"-01-02-00000.nc" "tmp/clm_"$year"0101.nc"
  ncks -v SOILWATER_10CM,TSOI_10CM,topo -d time,0,363 "../data/b40.20th.track1.1deg.012.clm2.h1."$year"-01-02-00000.nc" "tmp/clm_"$year"0102.nc"
  ncrcat "tmp/clm_"$year"0101.nc" "tmp/clm_"$year"0102.nc" "clm_"$year".nc"
  #CIC data
  ncrcat -v aice_d "../data/b40.20th.track1.1deg.012.cice.h1."$year* "cic_"$year".nc"
  ncrename -v TLAT,lat -v TLON,lon "cic_"$year".nc"
  ncatted -O -a coordinates,aice_d,o,c,"lon lat time" "cic_"$year".nc"
  #POP data
  ncrcat -v SST "../data/b40.20th.track1.1deg.012.pop.h.nday1."$year* "pop_"$year".nc"
  ncrename -v SST,tos -v TLAT,lat -v TLONG,lon "pop_"$year".nc"
  ncatted -O -a coordinates,tos,o,c,"lon lat time" "pop_"$year".nc"
  #CAM data
  ncrcat "../data/b40.20th.track1.1deg.012.cam2.h3."$yearb"-12"* "../data/b40.20th.track1.1deg.012.cam2.h3."$year* "tmp/cam_"$year".nc"
  cdo -selyear,$year "tmp/cam_"$year".nc" "cam_"$year".nc"
  ncks -A -v P0 "tmp/cam_"$year".nc" "cam_"$year".nc"

rm tmp/*.nc

done
