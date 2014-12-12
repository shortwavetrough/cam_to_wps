#!/bin/bash

cd $WORKDIR
months=(01 02 03 04 05 06 07 08 09 10 11 12)

for year in {1979..1989}; do
  for mon in {0..11}; do
    month=${months[$mon]}
    sed -e 's/st_yr/'${year}'/g' \
        -e 's/st_mo/'${month}'/g' \
        < CAM_netcdf_to_GRIB_tmpl.ncl > CAM_netcdf_to_GRIB.ncl
    ncl CAM_netcdf_to_GRIB.ncl &
    wait %1
  done
done

