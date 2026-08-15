
############################### from job.base.sh ###############################

        DATE1=`date "+%Y%m%d-%H:%M:%S"`
        OLD="old_${DATE1}"
        printf "\n date : ${DATE1}\n\n"


#===============================================================================
#  Step 1 commands:  get_file step  
#===============================================================================
if [ ${LOADL_STEP_NAME} == "get_file" ] || [ ${LOADL_STEP_NAME} == "XXX" ]; then
       
        printf " *****************************************************\n"
        printf " ****************** Step GET *************************\n"
        printf " *****************************************************\n\n" 
        #  Careful if modification in these vars, modif in the 3 steps of the job job.base.sh
        export JOBDIR="${JOBDIR_ROOT}/${ROOT_NAME_2}"
        export EXEDIR="${EXEDIR_ROOT}/${ROOT_NAME_2}"

        # EXEDIR: Execution directory
        [ -d ${EXEDIR} ] && mv ${EXEDIR} ${EXEDIR}_${OLD}
        mkdir -p ${EXEDIR} #&& mkdir ${EXEDIR}/ls_l

        # JOBDIR: job directory
        [ -d ${JOBDIR} ] && mv ${JOBDIR} ${JOBDIR}_${OLD}
        mkdir -p ${JOBDIR}

        printf " ************* Backup run files *****************\n\n"
        cpfile ${JOBDIR_ROOT}//myenv_mypath.sh ${JOBDIR}
        cpfile ${JOBDIR_ROOT}/mynamelist.sh ${JOBDIR}
        cpfile ${JOBDIR_ROOT}/myjob.sh ${JOBDIR}

        # some printings
        . ${SCRIPTDIR}/common_printing.sh

        cd ${EXEDIR}

       # Copy job in EXEDIR
       cp ${JOBDIR_ROOT}/${jobname} ./

#        printf "\n ************* SCRIPT files *****************\n\n"
#        for file in ${SCRIPTDIR}/*.sh ; do cpfile ${file} . ; done; chmod 755 *.sh

        printf "\n ************* Executable files *****************\n"
            if [ ${USE_OCE}  -eq 1 ]; then
                printf "\n  ***** OCEAN model *****\n\n" 
                if [ ${ONLINE_COMP} -eq 1 ]; then
                    . ${SCRIPTDIR}/oce_compile.sh 
                else
                    cpfile ${OCE_EXE_DIR}/croco.${RUNtype} crocox
                fi
                printf "\n   --- OCEAN model git version ---\n"
                . ${SCRIPTDIR}/getversion.sh ${OCE}
                printf "\n   -------------------------------\n"
            fi
	    if [ ${USE_ATM}  -eq 1 ]; then
                printf "\n  ***** ATMOSPHERIC model *****\n\n" 
                cpfile ${ATM_EXE_DIR}/wrf.exe wrfexe 
                printf "\n   --- ATMOSPHERIC model git version ---\n" 
                . ${SCRIPTDIR}/getversion.sh ${ATM}
                printf "\n   -------------------------------\n"
            fi
	    if [ ${USE_WAV}  -eq 1 ]; then
                printf "\n  ***** WAVE model *****\n\n" 
                printf "cp ${WAV_EXE_DIR}/ww3_* ."
                cp ${WAV_EXE_DIR}/ww3_* . 
                mv ww3_shel wwatch 
                printf "\n   --- WAVE model git version ---\n" 
                . ${SCRIPTDIR}/getversion.sh ${WAV}
                printf "\n   -------------------------------\n"
            fi

	    [ ${USE_XIOS} -ge 1 ] && { printf "\n  ***** XIOS server *****\n\n" ; cpfile ${XIOS_EXE_DIR}/xios_server.exe . ; }

            if [ ${USE_TOY}  -ge 1 ]; then
                printf "\n  ***** TOY model *****\n\n" 
		for k in `seq 0 $(( ${nbtoy} - 1))`; do
                    cpfile ${TOY_EXE_DIR}/toy_${toytype[$k]} toy${toytype[$k]}
		done
	    fi

        printf "\n *************** INPUT files ********************\n"

# if xios
        if [ ${USE_XIOS} -ge 1 ]; then
            printf "  ***** XIOS files *****\n\n"
	    [[ ${ONLINE_XML} == "TRUE" ]] && { . ${SCRIPTDIR}/xios_process.sh ;}
            cp ${XIOS_NAM_DIR}/*.xml . 
        fi

# ocean/atmosphere input files (configuration, forcing, obc, levitus/restart...)
        DATE_END_JOBm1=$( makedate $MONTH_BEGIN_JOB $(( $DAY_BEGIN_JOB - 1 )) $YEAR_BEGIN_JOB )
        RESTDIR_IN=${RESTDIR_ROOT}/${DATE_END_JOBm1}
        if [ ${USE_OCE} -eq 1 ]; then
            printf "\n  ***** get OCEAN CONFIGURATION, OBC, BLK... files *****\n\n"
            . ${SCRIPTDIR}/oce_getfile.sh
            printf "\n  ***** get OCEAN RESTART files *****\n\n" #|tee ls_l/oce_getrst.txt
            #printf "    see listing in ${EXEDIR}/ls_l/oce_getrst.txt"
            . ${SCRIPTDIR}/oce_getrst.sh #>> ls_l/oce_getrst.txt
            printf "\n  ***** prepare OCEAN namelist from *base file *****\n\n"
            . ${SCRIPTDIR}/oce_nam.sh
        fi
        if [ ${USE_ATM} -eq 1 ]; then
            printf "\n  ***** get ATMOSPHERE CONFIGURATION, INPUT, BDY... files *****\n\n" 
            . ${SCRIPTDIR}/atm_getfile.sh 
            printf "\n  ***** get ATMOSPHERE RESTART files *****\n\n" #|tee ls_l/atm_getrst.txt
            #printf "    see listing in ${EXEDIR}/ls_l/atm_getrst.txt"
            . ${SCRIPTDIR}/atm_getrst.sh #>> ls_l/atm_getrst.txt
            printf "\n  ***** prepare ATM namelist from $atmnamelist file *****\n\n"
            . ${SCRIPTDIR}/atm_nam.sh
        fi
        if [ ${USE_WAV} -eq 1 ]; then
            printf "\n  ***** get WAVE CONFIGURATION, inp, forcing ... files *****\n\n" 
            . ${SCRIPTDIR}/wav_getfile.sh 
            printf "\n  ***** prepare WAVE namelists *.inp file *****\n\n"
            . ${SCRIPTDIR}/wav_nam.sh
            printf "\n  ***** get WAVE RESTART files *****\n\n" #|tee ls_l/wav_getrst.txt
            #printf "    see listing in ${EXEDIR}/ls_l/wav_getrst.txt"
            . ${SCRIPTDIR}/wav_getrst.sh #>> ls_l/wav_getrst.txt
        fi
        echo $USE_CPL
        if [ ${USE_CPL} -ge 1 ]; then
            printf "\n  ***** get OASIS RESTART files *****\n\n" #|tee ls_l/cpl_getrst.txt
            #printf "    see listing in ${EXEDIR}/ls_l/cpl_getrst.txt"
            . ${SCRIPTDIR}/cpl_getrst.sh #>> ls_l/cpl_getrst.txt
            . ${SCRIPTDIR}/cpl_getfile.sh
            printf "\n  ***** prepare OASIS namelist from $namcouplename file *****\n\n"
            . ${SCRIPTDIR}/cpl_nam.sh
        fi
        if [ ${USE_TOY} -ge 1 ]; then
            printf "\n  ***** get TOY CONFIGURATION files *****\n\n"
            . ${SCRIPTDIR}/toy_getfile.sh
            printf "\n  ***** prepare TOY namelist from *base file *****\n\n"
            . ${SCRIPTDIR}/toy_nam.sh
        fi

        printf "\n date : `date "+%Y%m%d-%H:%M:%S"`\n"

fi # Step1
#===============================================================================
#  Step 2 commands:  run_model step (parallel)
#===============================================================================
if [ ${LOADL_STEP_NAME} == "run_model" ] || [ ${LOADL_STEP_NAME} == "XXX" ] ; then
       
        printf " ****************************************************\n"
        printf " ****************** Step RUN ************************\n"
        printf " ****************************************************\n\n"
 
        export JOBDIR="${JOBDIR_ROOT}/${ROOT_NAME_2}"
        export EXEDIR="${EXEDIR_ROOT}/${ROOT_NAME_2}"

cd ${EXEDIR} 

#       ls -l > ls_l/ls_pre_exe.txt
        #printf "\n see ls -l in ${EXEDIR}/ls_l/ls_pre_exe.txt\n"
        printf "\n date : `date "+%Y%m%d-%H:%M:%S"`\n"
        printf "\n---------------  start   ---------------\n"
#      
        printf "  make app.conf file for Multiple Program Multiple Data \n"
        [ -f ${SCRIPTDIR}/MACHINE/${MACHINE}/launch_${MACHINE}.sh ] && { . ${SCRIPTDIR}/MACHINE/${MACHINE}/launch_${MACHINE}.sh ; } || { printf "\n Please create a file launch_${MACHINE}.sh \n" ; exit 0 ; }

	printf "  launch run: $MPI_LAUNCH_CMD ${MPI_ext} app.conf "

##### RUN ######
time $MPI_LAUNCH_CMD ${MPI_ext} app.conf >& out_run.txt 
[ "$?" -gt "0" ] && { printf "ERROR during the RUN.\n Please check le log.files in ${EXEDIR} (out_run.txt, croco.log,...)"; exit 1 ; }
################

        printf "\n\n PWD: `pwd`\n"

        [ -f time.step ] && printf "\n  time.step: `cat time.step` \n"

        printf "\n---------------   end    ---------------\n"
        printf "\n date : `date "+%Y%m%d-%H:%M:%S"`\n"
        #printf "\n see ls -l in ${EXEDIR}/ls_l/ls_post_exe.txt\n"
        #ls -l > ls_l/ls_post_exe.txt
 
fi # Step2
 

#===============================================================================
#  Step 3 commands:  put_file step 
#===============================================================================
if [ ${LOADL_STEP_NAME} == "put_file" ] || [ "${LOADL_STEP_NAME}" == "XXX" ] ; then

        printf " ****************************************************\n"
        printf " ****************** Step PUT ************************\n"
        printf " ****************************************************\n\n"
     
        export JOBDIR="${JOBDIR_ROOT}/${ROOT_NAME_2}"
        export EXEDIR="${EXEDIR_ROOT}/${ROOT_NAME_2}"
        export OUTPUTDIR="${OUTPUTDIR_ROOT}/${ROOT_NAME_2}"
        export RESTDIR_OUT="${RESTDIR_ROOT}/${ROOT_NAME_3}"
        
# OUTPUTDIR: output directory
        ${MACHINE_STOCKAGE} ls ${OUTPUTDIR}  >  /dev/null  2>&1
        [ "$?" -eq "0" ] && ${MACHINE_STOCKAGE} mv ${OUTPUTDIR} ${OUTPUTDIR}_${OLD}
        ${MACHINE_STOCKAGE} mkdir -p ${OUTPUTDIR}

# RESTDIR_OUT: restart directory
        ${MACHINE_STOCKAGE} ls ${RESTDIR_OUT}  >  /dev/null  2>&1
        [ "$?" -eq "0" ] && ${MACHINE_STOCKAGE} mv ${RESTDIR_OUT} ${RESTDIR_OUT}_${OLD}
        ${MACHINE_STOCKAGE} mkdir -p ${RESTDIR_OUT}

cd ${EXEDIR} 
        printf "\n date : `date "+%Y%m%d-%H:%M:%S"`\n"

        if [ ${USE_OCE} -eq 1 ]; then
            printf "\n  ***** put OCEAN OUTPUT/RESTART files *****\n\n" #|tee ls_l/oce_putfile.txt
            #printf "    see listing in ${EXEDIR}/ls_l/oce_putfile.txt"
            . ${SCRIPTDIR}/oce_putfile.sh #>> ls_l/oce_putfile.txt
        fi
        if [ ${USE_ATM} -eq 1 ]; then
            printf "\n  ***** put ATMOSPHERE OUTPUT/RESTART files *****\n\n" #|tee ls_l/atm_putfile.txt
            #printf "    see listing in ${EXEDIR}/ls_l/atm_putfile.txt"
            . ${SCRIPTDIR}/atm_putfile.sh #>> ls_l/atm_putfile.txt
        fi
        if [ ${USE_WAV} -eq 1 ]; then
            printf "\n  ***** put WAVE OUTPUT/RESTART files *****\n\n" #|tee ls_l/wav_putfile.txt
            #printf "    see listing in ${EXEDIR}/ls_l/wav_putfile.txt"
            . ${SCRIPTDIR}/wav_putfile.sh #>> ls_l/wav_putfile.txt
        fi
        if [ ${USE_CPL} -ge 1 ]; then
            printf "\n  ***** put OASIS OUTPUT/RESTART files *****\n\n" #|tee ls_l/cpl_putfile.txt
            #printf "    see listing in ${EXEDIR}/ls_l/cpl_putfile.txt"
            . ${SCRIPTDIR}/cpl_putfile.sh #>> ls_l/cpl_putfile.txt
        fi 

        chmod -R ${permission} ${OUTPUTDIR} ${RESTDIR_OUT}
#-------------------------------------------------------------------------------
#  save output control ascii files in jobs directory
#-------------------------------------------------------------------------------

        printf "\n ***** save ascii job files in ${JOBDIR} *****\n\n" 
# if ocean
        FILES_OCE="layout.dat ocean.output* namelist out_run.txt *time.step solver.stat* croco.log"        
        [ ${USE_OCE} -eq 1 ] && {  for file in ${FILES_OCE}; do cpfile2 ${file} ${JOBDIR}; done; echo ""; }
# if atmosphere
        FILES_ATM="namelist.input"
        [ ${USE_ATM} -eq 1 ] && {  for file in ${FILES_ATM}; do cpfile2 ${file} ${JOBDIR}; done; echo ""; }
# if wave
        FILES_WAV="log.ww3 strt.out ounf.out prnc.*.out grid.out"
        [ ${USE_WAV} -eq 1 ] && {  for file in ${FILES_WAV}; do cpfile2 ${file} ${JOBDIR}; done; echo ""; }
# if coupler
        FILES_CPL="nout.000000 wrfexe.timers* crocox.timers* wwatch.timers*"
        [ ${USE_CPL} -ge 1 ] && {  for file in ${FILES_CPL}; do cpfile2 ${file} ${JOBDIR}; done; echo ""; }
# if agrif
        FILES_AGRIFZ="AGRIF_FixedGrids.in"
        [ ${AGRIFZ} -eq 1 ] && {  for file in ${FILES_AGRIFZ}; do cpfile2 ${file} ${JOBDIR}; done; echo ""; }
# job
        FILES_JOB="${jobname}"
# as for the moment steps are not in separated jobs, the log is not available as the job is not finised yet => commented
#        if [ ${MACHINE} == "IRENE" ]; then
#            FILES_JOB="${FILES_JOB} ${ROOT_NAME_1}*.o ${ROOT_NAME_1}*.e"
#        elif [ ${MACHINE} == "JEANZAY" ] || [ ${MACHINE} == "LEFTRARU" ]; then
#            FILES_JOB="${FILES_JOB} ${ROOT_NAME_1}.out"
#        elif [ ${MACHINE} == "DATARMOR" ] || [ ${MACHINE} == "WCHPC" ]; then
#            FILES_JOB="${FILES_JOB} ${ROOT_NAME_1}.o* ${ROOT_NAME_1}.e*"
#        else
#            FILES_JOB="${FILES_JOB} ${ROOT_NAME_1}.jobid_*.txt ${ROOT_NAME_1}.o*"
#        fi
        cd ${JOBDIR_ROOT}; for file in ${FILES_JOB}; do mvfile2 ${file} ${JOBDIR}; done;  cd -; echo "";

#  mv output files in output
        printf "\n date : `date "+%Y%m%d-%H:%M:%S"`\n"

#-------------------------------------------------------------------------------
#  NEXT job!
#-------------------------------------------------------------------------------

if [ "${MODE_TEST}" == "" ] ; then      #  en production
        printf "\n *************  run the next job  *****************\n\n"

        sed -e "s/YEAR_BEGIN_JOB=.*/YEAR_BEGIN_JOB=${YEAR_BEGIN_JOBp1}/" \
            -e "s/MONTH_BEGIN_JOB=.*/MONTH_BEGIN_JOB=${MONTH_BEGIN_JOBp1}/" \
            -e "s/DAY_BEGIN_JOB=.*/DAY_BEGIN_JOB=${DAY_BEGIN_JOBp1}/" \
            -e "s/NBJOB=.*/NBJOB=$(( ${NBJOB} - 1 ))/" \
            -e "s/RESTART_FLAG=.*/RESTART_FLAG=\"TRUE\"/" \
            -e "s/CHAINED_JOB=.*/CHAINED_JOB=\"FALSE\"/" \
                ${JOBDIR_ROOT}/myjob.sh > myjobtmp.sh
        mvfile2 ${JOBDIR_ROOT}/myjob.sh ${JOBDIR_ROOT}/myjob.sh.bck
        mvfile myjobtmp.sh ${JOBDIR_ROOT}/myjob.sh
        chmod 755   ${JOBDIR_ROOT}/myjob.sh

	if [ ${DATE_END_JOB} -ne ${DATE_END_EXP} ] ; then
            cd ${JOBDIR_ROOT}
            chmod 755 submitjob.sh
            ./submitjob.sh
        else
	    printf "\n ************* run finished at date : ${DATE_END_EXP} \n\n\n\n"
	    exit 0
	fi

else # en test
        printf "\n\n\n\n  MODE_TEST=${MODE_TEST}  Test mode and non production => No job chaining.\n\n\n\n"
fi

fi # Step3

