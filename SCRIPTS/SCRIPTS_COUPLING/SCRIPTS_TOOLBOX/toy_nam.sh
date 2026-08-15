. ${SCRIPTDIR}/caltools.sh
##
## -------------------------- #
## - Copy TOY input files -   #
## -------------------------- #

#-------
## Number of time step per day
#-------

echo "fill TOY model namelist"

for k in `seq 0 $(( ${nbtoy} - 1 ))` ; do
    if [ -f ${TOY_NAM_DIR}/${toynamelist[$k]} ] ; then
        #-------
	## Number of time step per day
	#-------
        TOY_NDT_DAY=$(( 86400 / ${DT_TOY[$k]} ))
	TOY_NTIMES=$(( ( ${JDAY_END_JOB} - ${JDAY_BEGIN_JOB} + 1 ) * ${TOY_NDT_DAY}     ))
        #
        sed -e "s/<toydt>/${DT_TOY[$k]}/g" -e "s/<toytimes>/${TOY_NTIMES}/g"   \
        ${TOY_NAM_DIR}/${toynamelist[$k]} > ./TOYNAMELIST.nam.${toytype[$k]}
    else
        printf "ERROR: could not find ${TOY_NAM_DIR}/${toynamelist[$k]}. Exit \n"
        exit 1
    fi
done
