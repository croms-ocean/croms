#!/bin/bash
#-------------------------------------------------------------------------------
#                                                          Configuration files
#-------------------------------------------------------------------------------

# Get grid file
for nn in $( seq 0 ${AGRIFZ} ); do
    if [ ${nn} -gt 0 ];    then
        agrif_ext=".${nn}"
    else
        agrif_ext=""
    fi
    printf "Get ocean grid file(s)\n"
    ${io_getfile} ${OCE_FILES_DIR}/croco_grd.nc${agrif_ext}                croco_grd.nc${agrif_ext}
    if [ ${nn} -gt 0 ];    then
        printf "Nesting: create AGRIF_FixedGrids.in file from infos contained in grid files"
        nest_pos=$( ncdump -v grd_pos croco_grd.nc${agrif_ext} | grep 'grd_pos =' | cut -d '=' -f 2 | cut -d ';' -f 1 )
        nest_coef=$( ncdump -v refine_coef croco_grd.nc${agrif_ext} | grep 'refine_coef =' | cut -d '=' -f 2 | cut -d ';' -f 1 )
        echo " 1" >> AGRIF_FixedGrids.in
        echo "$nest_pos $nest_coef $nest_coef $nest_coef $nest_coef" >> AGRIF_FixedGrids.in
        if [ ${nn} -eq ${AGRIFZ} ]; then
            echo " 0" >> AGRIF_FixedGrids.in
            sed -e "s/,/ /g" AGRIF_FixedGrids.in > AGRIF_FixedGrids.in.tmp
            mv AGRIF_FixedGrids.in.tmp AGRIF_FixedGrids.in
        fi
    fi
done

# Get boundary foring
printf "Get ocean boundary file \n"
. ${SCRIPTDIR}/oce_getbdy.sh

# Get surface forcing if needed
[ ${surfrc_flag} == "TRUE" ] && { printf "Get ocean surface forcing file(s) \n"; . ${SCRIPTDIR}/oce_getfrc.sh ; }

# Get tide forcing if any
for nn in $( seq 0 ${AGRIFZ} ); do
    if [ ${nn} -gt 0 ];    then
        agrif_ext=".${nn}"
    else
        agrif_ext=""
    fi
    [ ${tide_flag} == "TRUE" ] && { printf "Get tide forcing file(s)\n"; ${io_getfile} ${OCE_FILES_DIR}/croco_frc.nc${agrif_ext} ./ ; }
    [ ${river_flag} == "TRUE" ] && { printf "Get river forcing file(s)\n"; ${io_getfile} ${OCE_FILES_DIR}/croco_runoff.nc${agrif_ext} ./ ; }
done

          
	       














