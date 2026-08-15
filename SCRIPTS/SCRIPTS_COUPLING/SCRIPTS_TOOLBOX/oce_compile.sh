#!/bin/bash
if [[ ${RESTART_FLAG} == "FALSE" ]] ; then # || [[ ! -f "${OCE_EXE_DIR}/croco.${RUNtype}" ]]; then

#-------------------------------------------------------------------------------
#   Get files
#-------------------------------------------------------------------------------

    printf "   CROCO online compilation is on \n" 

#-------------------------------------------------------------------------------
#   sed on files
#-------------------------------------------------------------------------------
    cd ${JOBDIR_ROOT}
    cpfile ${OCE_EXE_DIR}/jobcomp .
    cpfile ${OCE_EXE_DIR}/cppdefs_dev.h .
    cpfile ${OCE_EXE_DIR}/param.h.base .
    cpfile ${OCE_EXE_DIR}/cppdefs.h.base .

    #   Jobcomp
    #----------
    printf "      Preparing jobcomp \n"
    sed -e "s|SOURCE=.*|SOURCE=${OCE} |g" \
        -e "s|SCRDIR=./Compile|SCRDIR=${JOBDIR_ROOT}/Compile |g" \
        -e "s|RUNDIR=.*|RUNDIR=${JOBDIR_ROOT} |g" \
        -e "s|FC=gfortran|FC=${FC}|g" \
        -e "s|MPIF90=.*|MPIF90=${MPIF90}|g" \
        -e "s|-O3|-O2|g" \
        ./jobcomp > tmp$$
    mv tmp$$ jobcomp
    #
    line=$(grep -n -m1 'PRISM_ROOT_DIR=' jobcomp | cut -d: -f1)
    sed "${line} s|PRISM_ROOT_DIR=.*|PRISM_ROOT_DIR=${CPL}|" ./jobcomp > tmp$$
    mv tmp$$ jobcomp
    #
    line=$(grep -n -m1 'XIOS_ROOT_DIR=' jobcomp | cut -d: -f1)
    sed "${line} s|XIOS_ROOT_DIR=.*|XIOS_ROOT_DIR=${XIOS}|" ./jobcomp > tmp$$
    mv tmp$$ jobcomp

    printf "      Reading grid size in ${OCE_FILES_DIR}/croco_grd.nc \n"
    cur_Y=$( echo $DATE_BEGIN_JOB | cut -c 1-4 )
    cur_M=$( echo $DATE_BEGIN_JOB | cut -c 5-6 )
    dimx=$( ncdump -h  ${OCE_FILES_DIR}/croco_grd.nc  | grep "xi_rho =" | cut -d ' ' -f 3)
    dimy=$( ncdump -h  ${OCE_FILES_DIR}/croco_grd.nc | grep "eta_rho =" | cut -d ' ' -f 3)
    dimz=$( ncdump -h  ${OCE_FILES_DIR}/croco_${ini_ext}_Y${cur_Y}M${cur_M}.nc | grep "s_rho =" | head -n 1 | cut -d ' ' -f 3)
    printf "          Grid size is (in Lx X Ly X Nz ) : ${dimx}X${dimy}X${dimz}\n"

    #   param.h
    #-----------
    printf "      Preparing param.h \n"
    # replace else conf settings in param.h
    sed -e "s/(\s*LLm0=xx,\s*MMm0=xx,\s*N=xx)/(LLm0=$(( ${dimx} - 2 )), MMm0=$(( ${dimy} - 2 )), N=${dimz})/g" \
        param.h.base > tmp$$
    mv tmp$$ param.h
    # update necessary things
    sed -e "s/NP_XI *= *[0-9]* *,/NP_XI=${NP_OCEX},/g" \
        -e "s/NP_ETA *= *[0-9]* *,/NP_ETA=${NP_OCEY},/g" \
        param.h > tmp$$ 
    mv tmp$$ param.h
    if [[ ${MPI_NOLAND} == "TRUE" ]]; then
      sed -e "s|NNODES=NP_XI\*NP_ETA|NNODES=${NP_OCE}|g" \
          param.h > tmp$
      mv tmp$$ param.h
    fi

    #   cppdefs.h
    #------------
    printf "      Preparing cppdefs.h \n"
    sed -e "s/# define BENGUELA_LR/# define ${CEXPER}/g" \
        -e "s/# undef  MPI/# define  MPI/g" \
        ./cppdefs.h.base > tmp$$
    mv tmp$$ cppdefs.h

    if [[ ${MPI_NOLAND} == "TRUE" ]]; then
      sed -e "s/# *undef  MPI_NOLAND/# define MPI_NOLAND/g" cppdefs.h > tmp$$
      mv tmp$$ cppdefs.h
    fi

    if [[ ${EXACT_RESTART} == "TRUE" ]]; then
      sed -e "s/# *undef *EXACT_RESTART/# define EXACT_RESTART/g" cppdefs.h > tmp$$
      mv tmp$$ cppdefs.h
    fi

    sed -e "s/# undef  LOGFILE/# define  LOGFILE/g" cppdefs.h > tmp$$
    mv tmp$$ cppdefs.h

    if [ $USE_CPL -ge 1 ]; then
        if [ $USE_ATM -eq 1 ] || [ $USE_TOYATM -eq 1 ]; then 
            sed -e "s/#  *undef  *OA_COUPLING/# define OA_COUPLING/g" cppdefs.h > tmp$$
	    mv tmp$$ cppdefs.h
            printf "           Coupling with ATM \n"
	    if [ ${CPL_PATM} == "TRUE" ]; then
	        printf "             Coupling PATM is defined \n"
                sed -e "s/# *undef *READ_PATM/#  define READ_PATM/g" \
                    -e "s/# *undef *OBC_PATM/#  define OBC_PATM/g" \
                    cppdefs.h > tmp$$
	    else
		printf "             Coupling PATM is NOT defined \n"
		printf "cppdefs.h is : "
                sed -e "s/# *define *READ_PATM/#  undef  READ_PATM/g" \
                    -e "s/# *define *OBC_PATM/#  undef  OBC_PATM/g" \
                    cppdefs.h > tmp$$
   	    fi
            mv tmp$$ cppdefs.h
	else
            sed -e "s/#  *define  *OA_COUPLING/# undef OA_COUPLING/g" cppdefs.h > tmp$$
	    mv tmp$$ cppdefs.h
	fi
        if [ $USE_WAV -eq 1 ] || [ $USE_TOYWAV -eq 1 ]; then
            printf "           Coupling with WAV \n"
            #mycase=`echo $RUNtype | tail -c 5`
            [[ $RUNtype =~ .*full.* ]] && mycase=full
            if [ $mycase == 'full' ]; then
                if [ ${OW_COUPLING_FULL} == "TRUE" ]; then 
                    sed -e "s/#  *undef  *OW_COUPLING/# define OW_COUPLING/g" \
                        -e "s/# *undef *MRL_WCI/# define MRL_WCI/g" \
                        cppdefs.h > tmp$$
                    printf "             OW_COUPLING_FULL option is activated \n"
                else
                    echo "ERROR... RUNtype contains 'full' but OW_COUPLING_FULL is not activated in mynamelist.sh. Inconsistent options. Exit"
                    exit 1
                fi
            else 
                if [[ ${OW_COUPLING_FULL} == "TRUE" ]]; then 
                    echo "ERROR... RUNtype does not contain 'full' but OW_COUPLING_FULL is activated in mynamelist.sh. Inconsistent options. Exit"
                    exit 1
                else
                    sed -e "s/#  *undef  *OW_COUPLING/# define OW_COUPLING/g" \
                        -e "s/# *define *OW_COUPLING_FULL/# undef OW_COUPLING_FULL/g" \
                        -e "s/# *undef *MRL_WCI/# define MRL_WCI/g" \
                        cppdefs.h > tmp$$
                fi
            fi
            mv tmp$$ cppdefs.h
            if [[ ${WAVE_SMFLUX} == "TRUE" ]]; then
                sed -e "s/#  *undef  *WAVE_SMFLUX/# define WAVE_SMFLUX/g" \
                    cppdefs.h > tmp$$
                printf "           WAVE_SMFLUX option is activated \n"
	        mv tmp$$ cppdefs.h
            fi
        else
            sed -e "s/#  *define  *OW_COUPLING/# undef OW_COUPLING/g" \
                -e "s/# *define *MRL_WCI/# undef  MRL_WCI/g" \
                cppdefs.h > tmp$$
	    mv tmp$$ cppdefs.h
        fi
            
    fi

    if [[ ${surfrc_flag} == "TRUE" && ${frc_ext} != *'frc'* ]]; then
	sed -e "s/#  *undef  *BULK_FLUX/# define BULK_FLUX/g" cppdefs.h > tmp$$
        mv tmp$$ cppdefs.h
        if [ ${interponline} == 1 ]; then
            sed -e "s/#  *undef  *ONLINE/# define ONLINE/g" cppdefs.h > tmp$$
            mv tmp$$ cppdefs.h
            if [[ ${frc_ext} == *'AROME'* || ${frc_ext} == *'ARPEGE'* ]]; then
                sed -e "s/#  *undef  *AROME/#   define AROME/g" \
                    cppdefs.h > tmp$$
            else
                sed -e "s/#  *define  *AROME/#   undef  AROME/g" \
                    cppdefs.h > tmp$$
            fi
            mv tmp$$ cppdefs.h
            if [[ ${frc_ext} == *'ERA_ECMWF'* ]]; then
                sed -e "s/#  *undef  *ERA_ECMWF/#   define  ERA_ECMWF/g" \
                    cppdefs.h > tmp$$
            else
                sed -e "s/#  *define  *ERA_ECMWF/#   undef  ERA_ECMWF/g" \
                    cppdefs.h > tmp$$
            fi
            printf "           Online bulk activated with ${frc_ext}\n"
            mv tmp$$ cppdefs.h
        elif [ ${interponline} == 0 ]; then
            sed -e "s/#  *define  *ONLINE/#  undef ONLINE/g" cppdefs.h > tmp$$
            mv tmp$$ cppdefs.h
            printf "           Bulk activated\n"
        fi
    elif [[ ${surfrc_flag} == "TRUE" && ${frc_ext} == *'frc'* ]]; then
        sed -e "s/#  *define  *BULK_FLUX/# undef BULK_FLUX/g"  cppdefs.h > tmp$$
    elif [[ ${surfrc_flag} == "FALSE" && $USE_ATM == 0 && $USE_TOYATM == 0 ]]; then
        echo "ERROR... surfrc_flag=FALSE and not coupling with atm is define. srfrc_flag should be TRUE. Exit"
        exit 1
    fi

    if [[ ${bdy_ext} == *'bry'* ]]; then
        sed -e "s/#  *define  *CLIMATOLOGY/# undef CLIMATOLOGY/g" \
	    -e "s/#  *undef *FRC_BRY/# define FRC_BRY/g" \
	cppdefs.h > tmp$$
        printf "           Lateral forcing is BRY\n"
        mv tmp$$ cppdefs.h
    else
        sed -e "s/#  *undef  *CLIMATOLOGY/# define CLIMATOLOGY/g" \
            -e "s/#  *define *FRC_BRY/# undef FRC_BRY/g" \
        cppdefs.h > tmp$$
        printf "           Lateral forcing is CLM\n"
    fi

    if [ ${tide_flag} == "TRUE" ]; then
	sed -e "s/#  *undef  *TIDES/# define TIDES/g" cppdefs.h > tmp$$
        mv tmp$$ cppdefs.h
        sed -e "s/#  *define  *TIDES_MAS/# undef TIDES_MAS/g" cppdefs.h > tmp$$
        mv tmp$$ cppdefs.h
        printf "           Tides are taken into account\n"
    else
        sed -e "s/#  *define  *TIDES/# undef TIDES/g" cppdefs.h > tmp$$
        mv tmp$$ cppdefs.h
    fi

    if [ ${river_flag} == "TRUE" ]; then
        sed -e "s/#  *undef  *PSOURCE/# define PSOURCE/g" \
            -e "s/#  *undef  *PSOURCE_NCFILE/# define PSOURCE_NCFILE/g" \
            -e "s/#  *undef *PSOURCE_NCFILE_TS/#  define PSOURCE_NCFILE_TS/g" \
            cppdefs.h > tmp$$
        mv tmp$$ cppdefs.h
        printf "           Rivers are taken into account\n"
    else
        sed -e "s/#  *define  *PSOURCE/# undef PSOURCE/g"\
            -e "s/#  *define  *PSOURCE_NCFILE/# undef PSOURCE_NCFILE/g" \
            -e "s/#  *define  *PSOURCE_NCFILE_TS/#  undef PSOURCE_NCFILE_TS/g" \
             cppdefs.h > tmp$$
        mv tmp$$ cppdefs.h
    fi

    if [ ${USE_XIOS_OCE} -eq 1 ]; then
	sed -e "s/#  *undef  *XIOS/# define XIOS/g" cppdefs.h > tmp$$
    	mv tmp$$ cppdefs.h
        printf "           Output will be handled by XIOS\n"
        if [ ${USE_XIOS_ATM} -eq 1 ]; then
            sed -e "s/#  *undef  *XIOS_ATM/# define XIOS_ATM/g" cppdefs_dev.h > tmp$$
            mv tmp$$ cppdefs_dev.h
            printf "           XIOS for ATM is also on\n"
        else
            sed -e "s/#  *define  *XIOS_ATM/# undef XIOS_ATM/g" cppdefs_dev.h > tmp$$
            mv tmp$$ cppdefs_dev.h
        fi
    else
        sed -e "s/#  *define  *XIOS/# undef XIOS/g" cppdefs.h > tmp$$
        mv tmp$$ cppdefs.h
    fi

    if [ $AGRIFZ -eq 0 ]; then
        sed -e "s/#  *define  *AGRIF/# undef AGRIF/g" \
            -e "s/#  *define   *AGRIF_2WAY/# undef AGRIF_2WAY/g" \
            cppdefs.h > tmp$$
        sed -e "s/MAKE  *\-j  *[1-9]/MAKE -j 8/g" jobcomp > tmp2$$
        
    else
        sed -e "s/#  *undef  *AGRIF/# define AGRIF/g" cppdefs.h > tmp$$
        mv tmp$$ cppdefs.h
        printf "           AGRIF is activated\n"
        if [[ ${AGRIF_2WAY} == "TRUE" ]]; then
            sed -e "s/#  *undef  *AGRIF_2WAY/# define AGRIF_2WAY/g" cppdefs.h > tmp$$
            printf "           AGRIF_2WAY is activated\n"
        else
            sed -e "s/#  *define  *AGRIF_2WAY/# undef AGRIF_2WAY/g" cppdefs.h > tmp$$
	fi
        sed -e "s/MAKE  *\-j  *[1-9]/MAKE -j 1/g" jobcomp > tmp2$$
    fi

    mv tmp$$ cppdefs.h
    mv tmp2$$ jobcomp

#-------------------------------------------------------------------------------
#   compile
#-------------------------------------------------------------------------------
    printf "     Launching jobcomp... \n"
    chmod 755 jobcomp
    time ./jobcomp >& log.compil
    if [ -f croco ]; then
        printf "\n Successful compilation, copying executable in ${EXEDIR} \n"
    else
        printf "\nERROR while compiling CROCO.\n Please check ${PWD}/log.compil"
        exit 1
    fi
    mv croco croco.${RUNtype}
    cp cppdefs.h cppdefs.h.${RUNtype}
    cp param.h param.h.${RUNtype}
    # save exe for next jobs
    cpfile croco.${RUNtype} ${EXEDIR}/crocox
    #    [[ ${USE_XIOS_OCE} == 1 && -d "ls -A ${XIOS_NAM_DIR}" ]] && { cp *.xml ${XIOS_NAM_DIR}/ ;}
    cd ${EXEDIR}
else
    printf "\n Copying executable in ${EXEDIR} \n" 
    cpfile ${JOBDIR_ROOT}/croco.${RUNtype} crocox
    #    [[ ${USE_XIOS_OCE} == 1 && -d "ls -A ${XIOS_NAM_DIR}" ]] && { cp *.xml ${XIOS_NAM_DIR}/ ;}
    
fi

