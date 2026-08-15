#-------------------------------------------------------------------------------
#                                                                      Average
#-------------------------------------------------------------------------------

year=$( printf "%04d"   ${YEAR_BEGIN_JOB} )
month=$( printf "%02d"  ${MONTH_BEGIN_JOB} )


# WW3 format to netcdf
 if [ -e out_grd.ww3 ]; then
  echo 'WW3 post-processing of gridded data:'
  echo "${SERIAL_LAUNCH_WAV}ww3_ounf &> ounf.out"
  ${SERIAL_LAUNCH_WAV}ww3_ounf &> ounf.out 

  #module load $ncomod
  #ncrcat -O ww3.??????.nc ww3.${year}${month}.nc #Concatenate all ww3*.nc file before exporting to OUTPUTDIR
  #ncrcat -O -F -d time,1,-2 ww3.${year}${month}.nc ww3.${year}${month}.nc
  #module unload $ncomod

  #${io_putfile} ww3.${year}${month}.nc ${OUTPUTDIR}/ww3_${DATE_BEGIN_JOB}_${DATE_END_JOB}.nc

 if [ "${flagout}" == "TRUE" ]; then
     echo "Keep WW3 output binary file out_grd.ww3"
     ${io_putfile} out_grd.ww3 ${OUTPUTDIR}/out_grd_${DATE_BEGIN_JOB}_${DATE_END_JOB}.ww3
 fi
 fi
 
 if [ -e out_pnt.ww3 ]; then
  echo 'WW3 post-processing of point output (spectra):'
  # WW3 format to netcdf
  echo "${SERIAL_LAUNCH_WAV}ww3_ounp &> ounp.out"
  ${SERIAL_LAUNCH_WAV}ww3_ounp &> ounp.out

  #module load $ncomod
  #ncrcat -O ww3.??????_spec.nc ww3.${year}${month}_spec.nc #Concatenate all ww3*.nc file before exporting to OUTPUTDIR
  #ncrcat -O -F -d time,1,-2 ww3.${year}${month}_spec.nc ww3.${year}${month}_spec.nc
  #module unload $ncomod

  #${io_putfile} ww3.${year}${month}_spec.nc ${OUTPUTDIR}/ww3_${DATE_BEGIN_JOB}_${DATE_END_JOB}_spec.nc

 if [ "${flagout}" == "TRUE" ]; then
     echo "Keep WW3 output binary file out_pnt.ww3"
     ${io_putfile} out_pnt.ww3 ${OUTPUTDIR}/out_pnt_${DATE_BEGIN_JOB}_${DATE_END_JOB}.ww3
 fi
 fi

mv ww3.*.nc ${OUTPUTDIR}/.

#-------------------------------------------------------------------------------
#                                                                      Restart
#-------------------------------------------------------------------------------

rstfile='mod_def.ww3 restart001.ww3'
for k in `seq 0 $(( ${lengthforc} - 1))` ; do
    [[ `echo ${forcww3[$k]} | head -c 4` == 'wind' ]] && forctype=wind
    [[ `echo ${forcww3[$k]} | head -c 3` == 'cur' ]] && forctype=current
    [[ `echo ${forcww3[$k]} | head -c 3` == 'lev' ]] && forctype=level
    rstfile="${rstfile} ${forctype}.ww3"
done

for file in ${rstfile} ; do 
    if [ ${file} = 'restart001.ww3' ]; then
        filesend='restart.ww3'
    else
        filesend=${file}
    fi
    ${io_putfile} ${file} ${RESTDIR_OUT}/${filesend}_${DATE_END_JOB}
done
 
