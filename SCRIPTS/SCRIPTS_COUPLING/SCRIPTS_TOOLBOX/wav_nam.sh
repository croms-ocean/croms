#!/bin/bash

echo 'Copy *.inp files'
cp ${WAV_NAM_DIR}/*.inp .

ms=$( printf "%02d"  ${MONTH_BEGIN_JOB} )
me=$( printf "%02d"  ${MONTH_END_JOB} )
ds=$( printf "%02d"  ${DAY_BEGIN_JOB} )
de=$( printf "%02d"  ${DAY_END_JOB} )

echo 'Fill ww3_grid.inp file'
sed -e "s/<wavdt>/${DT_WAV}/g" \
    -e "s/<wavdtPRO>/${DT_WW_PRO}/g"  -e "s/<wavdtREF>/${DT_WW_REF}/g"  -e "s/<wavdtSRC>/${DT_WW_SRC}/g"  \
    -e "s/<wavnx>/${wavnx}/g"   -e "s/<wavny>/${wavny}/g"  \
    -e "s/<hmin>/${hmin}/g" \
    -e "s/<CEXPER>/${CEXPER}/g" \
    ${WAV_NAM_DIR}/ww3_grid.inp.base > ./ww3_grid.inp

echo 'Fill ww3_shel.inp file'
if [ -f ${point_output_list} ] ; then
    sed 's/^/   /' ${point_output_list} > tmp_point_output_list
    echo "   0.0   0.0  'STOPSTRING'" | cat >> tmp_point_output_list
    #my_point_output_list="tmp_point_output_list"
fi
if [ ${wav_trck} -ne 0 ]; then 
    my_track_flag='   T'
else
    my_track_flag='$'
fi
sed -e "s/<yr1>/${YEAR_BEGIN_JOB}/g"  -e "s/<mo1>/${ms}/g"  -e "s/<dy1>/${ds}/g"  -e "s/<hr1>/00/g"  \
    -e "s/<yr2>/${YEAR_END_JOB}/g"  -e "s/<mo2>/${me}/g"  -e "s/<dy2>/${de}/g"  -e "s/<hr2>/24/g" \
    -e "s/<wav_int>/${wav_int}/g"  -e "s/<wav_rst>/$(( ${TOTAL_JOB_DUR} * 24 * 3600))/g" \
    -e "s/<wav_pnt>/${wav_pnt}/g"  -e "s/<wav_trck>/${wav_trck}/g" \
    -e "/$     <point_output_list>/r tmp_point_output_list"  \
    -e "s/$    <T>/$my_track_flag/g" \
    -e "s/<wavdt>/${DT_WAV}/g" \
    ${WAV_NAM_DIR}/ww3_shel.inp.base.${SHELL_EXT} > ./ww3_shel.inp

echo 'Fill ww3_ounf.inp file'
sed -e "s/<wav_int>/${wav_int}/g" \
    -e "s/<yr1>/${YEAR_BEGIN_JOB}/g"  -e "s/<mo1>/${ms}/g"  -e "s/<dy1>/${ds}/g"  -e "s/<hr1>/00/g" \
    ${WAV_NAM_DIR}/ww3_ounf.inp.base > ./ww3_ounf.inp

if [ ${wav_pnt} -ne 0 ]; then
    if [ -f ${point_output_list} ] ; then  
        echo 'Fill ww3_ounp.inp file'
        sed -e "s/<wav_pnt>/${wav_pnt}/g" \
            -e "s/<yr1>/${YEAR_BEGIN_JOB}/g"  -e "s/<mo1>/${ms}/g"  -e "s/<dy1>/${ds}/g"  -e "s/<hr1>/00/g" \
            ${WAV_NAM_DIR}/ww3_ounp.inp.base > ./ww3_ounp.inp
    else
         echo "ERROR: wav_pnt is not 0 and point_output_list=${point_output_list} file cannot be found. Exit"
         exit 1
    fi
fi

if [ ! -z $bouncin ]; then
    echo 'Fill the ww3_bounc.inp file'
    for bdyfile in ${WAV_FILES_DIR}/$bouncin* ; do
        echo $(basename $bdyfile) >> spec_bdy_list
    done
    sed -e "/'STOPSTRING'/{ h
            e cat spec_bdy_list
            g
            }" ${WAV_NAM_DIR}/ww3_bounc.inp > ww3_bounc.inp
fi

#if [ ${wav_trck} -ne 0 ]; then 
#    echo 'Fill ww3_trnc.inp file'
#    sed -e "s/<wav_trck>/${wav_trck}/g" \
#        -e "s/<yr1>/${YEAR_BEGIN_JOB}/g"  -e "s/<mo1>/${ms}/g"  -e "s/<dy1>/${ds}/g"  -e "s/<hr1>/00/g" \
#        ${WAV_NAM_DIR}/ww3_trnc.inp.base > ./ww3_trnc.inp
#fi

