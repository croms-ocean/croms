#!/bin/bash
#
########################################################
#  Define files and run parameters
########################################################
#
# Name used for the input files. For example croco_grd.nc
MODEL=croco

# Scratch directory where the model is run
SCRATCHDIR=`pwd`/SCRATCH

# Input directory where the croco_inter.in input file is located
INPUTDIR=`pwd`/CROCO_IN  # prod architecture
#INPUTDIR=`pwd`          # dev architecture

# AGRIF input file which defines the position of child grids
AGRIF_FILE=AGRIF_FixedGrids.in

# Directory where the croco input NetCDF files (croco_grd.nc, ...) are stored
MSSDIR=`pwd`/CROCO_FILES

# Directory where the croco output and restart NetCDF files (croco_his.nc, ...) are stored
MSSOUT=$SCRATCHDIR

# CROCO executable
CODFILE=croco

# number of processors for MPI run
NBPROCS=8

# command for running the mode : ./ for sequential job, mpirun -np NBPROCS for mpi run
# WARNING: for mpi run command, it is needed to add a space at the end!
RUNCMD='./'
#RUNCMD="mpirun -np $NBPROCS "
#RUNCMD="$MPI_LAUNCH "
#RUNCMD='srun '
#RUNCMD="mpiexec -np ${NBPROCS} ./"

#  Define environment variables for OPENMP
OMP_SCHEDULE=static
OMP_NUM_THREADS=2
OMP_DYNAMIC=false
OMP_NESTED=false
KMP_LIBRARY=throughput
KMP_STACKSIZE=2m
KMP_DUPLICATE_LIB_OK=TRUE

# Define which type of inputs are used
BULK_FILES=1
FORCING_FILES=0
CLIMATOLOGY_FILES=0
BOUNDARY_FILES=1
RUNOFF_FILES=0
TIDE_FILES=0
ONLINE_FILES=0
ONLINEFREQ=4
ONLINEPATH="../DATA/CFSR_Benguela_LR" # we recommend using an absolute path

# Define the suffix of your input files
# Atmospheric bulk file suffix (croco_blk_ATMOS_BULK_Y????M??.nc) - NOT USED IF ONLINE
ATMOS_BULK=ERA5
# Atmospheric forcing file suffix (croco_frc_ATMOS_FRC_Y????M??.nc) - NOT USED IF BULK or ONLINE
ATMOS_FRC=QSCAT
# Oceanic boundary and initial files suffix (croco_ini_OGCM_Y????M??.nc and croco_bry_OGCM_Y????M??.nc)
OGCM=SODA
# Runoff file sufix (croco_runoff_RUNOFF_Y????M??.nc)
RUNOFF=DAI
# Tide file suffix (croco_frc_TIDE_FRC.nc)
TIDE_FRC=tpxo7_croco

# Model time step [seconds]
DT=3600
# Number of barotropic time steps within one baroclinic time step [number], NDTFAST in croco.in
NFAST=60

# Number total of grid levels (1: No child grid)
NLEVEL=1
# AGRIF nesting refinement coefficient
AGRIF_REF=3

# Start and End year 
NY_START=2005
NY_END=2005
# Start and End month
NM_START=1
NM_END=3
# Set month format at 1 or 2 digits (for input and output files): "%01d" = 1 digit/ "%02d" = 2 digit  
MTH_FORMAT="%02d"
#  Time Schedule  -  TIME_SCHED=0 --> yearly files
#                    TIME_SCHED=1 --> monthly files
TIME_SCHED=1

# Number of year that are considered to be part of the spin-up (i.e. 365 days per year)
NY_SPIN=0

# Output frequency [days] - case 1 : No USE_CALENDAR - DEFAULT
#   average
ND_AVG=3
#   history (if = -1 set equal to NUMTIMES, the end of each month/year)
ND_HIS=-1
#   restart (if = -1 set equal to NUMTIMES, the end of each month/year)
ND_RST=-1

# Output frequency [hours] - case 2 : USE_CALENDAR - USED ONLY IF USE_CALENDAR
#
USE_CALENDAR=1
#
#  average (in hours)
NHAVG_UC=$((24))
#  history (in hours, if = -1 set equal t NUMTIMES*DT/3600, the end of each month/year)
NHHIS_UC=-1
#  restart (in hours, if = -1 set equal to NUMTIMES*DT/3600, the end of each month/year)
NHRST_UC=-1

#  Restart file - RSTFLAG=0 --> No Restart
#		  RSTFLAG=1 --> Restart
RSTFLAG=0     
#  Exact restart - EXACT_RST=0 --> Exact restart OFF
#                - EXACT_RST=1 --> Exact restart ON
EXACT_RST=0

#unalias cp
#unalias mv
#limit coredumpsize unlimited
CP=/bin/cp
MV=/bin/mv
LN=/bin/ln
#
########################################################
#  END USER CHANGE
########################################################
#
if [[ $TIME_SCHED == 0 ]]; then
  NM_START=1979
  NM_END=1979
fi
#
# netcdf file prefixes
#
GRDFILE=${MODEL}_grd
FRCFILE=${MODEL}_frc
BLKFILE=${MODEL}_blk
INIFILE=${MODEL}_ini
CLMFILE=${MODEL}_clm
BRYFILE=${MODEL}_bry
RNFFILE=${MODEL}_runoff
TIDEFILE=${MODEL}_frc
#
if [ ! -e $MSSOUT ] ; then
 mkdir $MSSOUT
fi
#
if [[ $RSTFLAG != 0 ]]; then
  NY=$NY_START
  NM=$NM_START
  if [[ $TIME_SCHED == 0 ]]; then
    NY=$((NY - 1))
    TIME=Y${NY}
  else
    NM=$((NM - 1))
    if [[ $NM == 0 ]]; then
      NM=12
      NY=$((NY - 1))
    fi
    TIME=Y${NY}M$( printf ${MTH_FORMAT} ${NM})
  fi
  RSTFILE=${MODEL}_rst_${TIME}
fi
#
if [[ $TIME_SCHED == 0 ]]; then
  TIME=Y${NY_START}
else
  TIME=Y${NY_START}M$( printf ${MTH_FORMAT} ${NM_START})
fi
#
# Get the code
#
if [ ! -e $SCRATCHDIR ] ; then
 mkdir $SCRATCHDIR
fi
cd $SCRATCHDIR
echo "Getting $CODFILE from $INPUTDIR"
$CP -f $INPUTDIR/$CODFILE $SCRATCHDIR
chmod u+x $CODFILE
if [[ $NLEVEL > 1 ]]; then
  echo "Getting $AGRIF_FILE from $MSSDIR"
  $CP -f $MSSDIR/$AGRIF_FILE $SCRATCHDIR
fi
#
# Get the netcdf files
#
LEVEL=0
while [ $LEVEL != $NLEVEL ]; do
  if [[ ${LEVEL} == 0 ]]; then
    ENDF=
  else
    ENDF=.${LEVEL}
  fi
  echo "Getting ${GRDFILE}.nc${ENDF} from $MSSDIR"
  $LN -sf $MSSDIR/${GRDFILE}.nc${ENDF} $SCRATCHDIR
  echo "Getting ${MODEL}_inter.in${ENDF} from $INPUTDIR"
  $CP -f $INPUTDIR/${MODEL}_inter.in${ENDF} $SCRATCHDIR
  if [[ $RSTFLAG == 0 ]]; then
    echo "Getting ${INIFILE}_${OGCM}_${TIME}.nc${ENDF} from $MSSDIR"
    $CP -f $MSSDIR/${INIFILE}_${OGCM}_${TIME}.nc${ENDF} $SCRATCHDIR
    $CP -f ${INIFILE}_${OGCM}_${TIME}.nc${ENDF} ${INIFILE}.nc${ENDF}
  else
    echo "Getting ${RSTFILE}.nc${ENDF} from $MSSOUT"
    $CP -f $MSSOUT/${RSTFILE}.nc${ENDF} $SCRATCHDIR
    $CP -f ${RSTFILE}.nc${ENDF} ${INIFILE}.nc${ENDF}
  fi
  LEVEL=$((LEVEL + 1))
done
###########################################################
#  Compute
###########################################################
#
NY_END=$((NY_END + 1))
NM_END=$((NM_END + 1))
NY=$NY_START
while [ $NY != $NY_END ]; do
  if [[ $NY == $NY_START ]]; then
    NM=$NM_START
  else 
     NM=1
  fi
   MY_YEAR=$NY
   MY_YEAR=$((MY_YEAR + 1))
  if [[ $MY_YEAR == $NY_END ]]; then
     MONTH_END=$NM_END
  else 
     MONTH_END=13
  fi
  if [[ $TIME_SCHED == 0 ]]; then
     MONTH_END=2
  fi
  while [ $NM != $MONTH_END ]; do
    if [[ $TIME_SCHED == 0 ]]; then
      TIME=Y${NY}
      echo "Computing YEAR $NY"
    else
	TIME=Y${NY}M$( printf ${MTH_FORMAT} ${NM})
	echo "Computing YEAR $NY MONTH $( printf ${MTH_FORMAT} ${NM})"
    fi
#
# Get forcing and clim for this time
#
    LEVEL=0
    while [ $LEVEL != $NLEVEL ]; do
      if [[ ${LEVEL} == 0 ]]; then
        ENDF=
      else
        ENDF=.${LEVEL}
      fi
      if [[ ${FORCING_FILES} == 1 ]]; then
        echo "Getting ${FRCFILE}_${ATMOS_FRC}_${TIME}.nc${ENDF} from $MSSDIR"
        $LN -sf $MSSDIR/${FRCFILE}_${ATMOS_FRC}_${TIME}.nc${ENDF} ${FRCFILE}.nc${ENDF}
      fi
      if [[ ${BULK_FILES} == 1 ]]; then
          echo "Getting ${BLKFILE}_${ATMOS_BULK}_${TIME}.nc${ENDF} from $MSSDIR"
          $LN -sf $MSSDIR/${BLKFILE}_${ATMOS_BULK}_${TIME}.nc${ENDF} ${BLKFILE}.nc${ENDF}
      fi
      if [[ ${RUNOFF_FILES} == 1 ]]; then
          echo "Getting ${RNFFILE}.nc${ENDF} from $MSSDIR"
          $LN -sf $MSSDIR/${RNFFILE}.nc${ENDF} ${RNFFILE}.nc${ENDF}
      fi
      if [[ ${TIDE_FILES} == 1 ]]; then
          echo "Getting ${TIDEFILE}_${TIDE_FRC}.nc${ENDF} from $MSSDIR"
	  $LN -sf $MSSDIR/${TIDEFILE}_${TIDE_FRC}.nc${ENDF} ${TIDEFILE}.nc${ENDF}
      fi

      LEVEL=$((LEVEL + 1))
    done
#
# No child climatology or boundary files
#
    if [[ ${CLIMATOLOGY_FILES} == 1 ]]; then
      echo "Getting ${CLMFILE}_${OGCM}_${TIME}.nc from $MSSDIR"
      $LN -sf $MSSDIR/${CLMFILE}_${OGCM}_${TIME}.nc ${CLMFILE}.nc
    fi
    if [[ ${BOUNDARY_FILES} == 1 ]]; then
      echo "Getting ${BRYFILE}_${OGCM}_${TIME}.nc from $MSSDIR"
      $LN -sf $MSSDIR/${BRYFILE}_${OGCM}_${TIME}.nc ${BRYFILE}.nc
    fi
#
# Set the number of time steps for each month 
# (30 or 31 days + 28 or 29 days for february)
#
    NUMTIMES=0
#
    if [[ ${NM} == 1 || ${NM} == 3 || ${NM} == 5 || ${NM} == 7 || ${NM} == 8 || ${NM} == 10 || ${NM} == 12 ]]; then
      NDAYS=31
    else
      NDAYS=30
      if [[ ${NM} == 2 ]]; then
        NDAYS=28
# February... check if it is a leap year

        B4=0
        B100=0
        B400=0

        B4=$((4 * ( NY / 4 )))
        B100=$((100 * ( NY / 100 )))
        B400=$((400 * ( NY / 400 )))

	
        if [[ $NY == $B4 && ((!($NY == $B100))||($NY == $B400)) ]]; then
#
          BSPIN=$(( NY - NY_START + 1 ))
          if [[ $BSPIN -gt $NY_SPIN ]]; then
	     echo Leap Year - $NY $B4 $B100 $B400
             NDAYS=29
          else
#.........   SPINUP!!!! In case of spinup I cant have leap years.
	     echo year $NY should be a Leap Year     
	     echo 'BUT : Spinup case: no leap year'
             NDAYS=28
          fi
#
        else
	  echo Not a Leap Year - $NY $B4 $B100 $B400
          NDAYS=28	  		  
        fi
      fi
    fi
    #
    # Put the number of time steps in the .in files
    #
    echo "YEAR = $NY MONTH = $NM DAYS = $NDAYS DT = $DT NTIMES = $NUMTIMES"
    NUMTIMES=$((NDAYS * 24 * 3600))
    NUMTIMES=$((NUMTIMES / DT))

    DT0=$DT
    LEVEL=0
    while [[ $LEVEL != $NLEVEL ]]; do
	if [[ ${LEVEL} == 0 ]]; then
            ENDF=
	else
            ENDF=.${LEVEL}
	    NUMTIMES=$((AGRIF_REF * NUMTIMES))
	    DT=$((DT / AGRIF_REF))
	fi
	NUMAVG=$((ND_AVG * 86400 / DT ))
	if [[ ${ND_HIS} -ne -1 ]]; then
	    NUMHIS=$((ND_HIS * 86400 / DT ))
	else
	    NUMHIS=$NUMTIMES
	fi
	if [[ ${ND_RST} -ne -1 ]]; then
	    NUMRST=$((ND_RST * 86400 / DT ))
	else
	    NUMRST=$NUMTIMES
	fi
	if [[ $USE_CALENDAR == 1 ]]; then
	    echo "USE_CALENDAR defined"
	    NHAVG=$((NHAVG_UC))
	    if [[ ${NHHIS_UC} -ne -1 ]]; then
		NHHIS=$((NHHIS_UC))
	    else
		NHHIS=$((NDAYS * 24))
	    fi
	    if [[ ${NHRST_UC} -ne -1 ]]; then
		NHRST=$((NHRST_UC))
	    else
		NHRST=$((NDAYS * 24))
	    fi
	fi
	if [[ $EXACT_RST == 1 ]]; then
	    echo "Exact restart defined"
	    if [[ $NY == $NY_START && $NM == $NM_START ]]; then
		NUMRECINI=1
		echo "set NUMRECINI = $NUMRECINI"
	    else
		NUMRECINI=2
		echo "set NUMRECINI = $NUMRECINI"
	    fi
	else  # no exact restart
	    echo "No exact restart"
	    NUMRECINI=1
	    echo "set NUMRECINI = $NUMRECINI"
	fi
	
	echo " "
	echo "Writing in ${MODEL}_inter.in${ENDF}"
	echo "USING DT       = $DT"
	echo "USING NFAST    = $NFAST"
	echo "USING NUMTIMES = $NUMTIMES"
	echo "USING NUMAVG   = $NUMAVG"
	echo "USING NUMHIS   = $NUMHIS"
	echo "USING NUMRST   = $NUMRST"
	echo "USING NUMRECINI = $NUMRECINI"
	echo "USING NYONLINE = $NY"
	echo "USING NMONLINE = $NM"
	echo "USING ENDYONLINE = $NY_END"
	echo "USING ENDMONLINE = $NM_END"
	echo "USING ONLINEFREQ = $ONLINEFREQ"
	echo "USING ONLINEPATH = $ONLINEPATH"
	
	if [ ! -f ${MODEL}_inter.in${ENDF} ]; then
	    echo "=="
	    echo "=> ERROR : miss the ${MODEL}_inter.in${ENDF} file"
	  echo "=="
	  exit 1
	fi
        sed \
          -e "s/NUMTIMES/${NUMTIMES}/" \
          -e "s/TIMESTEP/${DT}/" \
          -e "s/NFAST/${NFAST}/" \
          -e "s/\bNUMAVG\b/${NUMAVG}/" \
          -e "s/\bNUMHIS\b/${NUMHIS}/" \
          -e "s/\bNUMRST\b/${NUMRST}/" \
          -e "s/NUMRECINI/${NUMRECINI}/" \
          -e "s/NYONLINE/${NY}/" \
          -e "s/NMONLINE/${NM}/" \
          -e "s/ENDYONLINE/${NY_END}/" \
          -e "s/ENDMONLINE/${NM_END}/" \
          -e "s/ONLINEFREQ/${ONLINEFREQ}/" \
          -e "s|ONLINEPATH|${ONLINEPATH}|" \
          -e "s|<logfilename>|${MODEL}_${TIME}.out|" \
          < "${MODEL}_inter.in${ENDF}" > "${MODEL}_${TIME}_inter.in${ENDF}"

	if [[ $USE_CALENDAR == 1 ]]; then
	    if [[ ${NM} == 12 ]]; then
		NM_E_UC=1
		NY_E_UC=$((NY + 1))
	    else
		NM_E_UC=$((NM + 1))
		NY_E_UC=$NY
	    fi
	    echo "USING Ystart   = $NY"
	    echo "USING Mstart   = $(printf "%02d" $NM)"
	    echo "USING Yend     = $NY_E_UC"
	    echo "USING Mend     = $(printf "%02d" $NM_E_UC)"
	    echo "USING NHHIS    = $NHHIS"
	    echo "USING NHAVG    = $NHAVG"
	    echo "USING NHRST    = $NHRST"
	    sed -e "s/NHHIS/${NHHIS}/" \
		-e "s/NHAVG/${NHAVG}/" \
		-e "s/NHRST/${NHRST}/" \
	    	-e "s/Ystart/${NY}/"   \
		-e "s/Mstart/$(printf "%02d" $NM)/" \
		-e "s/Yend/${NY_E_UC}/" \
		-e "s/Mend/$(printf "%02d" $NM_E_UC)/" \
		< "${MODEL}_${TIME}_inter.in${ENDF}" > "${MODEL}_${TIME}_inter_UC.in${ENDF}"
	    mv "${MODEL}_${TIME}_inter_UC.in${ENDF}" "${MODEL}_${TIME}_inter.in${ENDF}"
	    fi
	#
	LEVEL=$((LEVEL + 1))
    done
    DT=$DT0
    #
    #  COMPUTE
    #
    echo " "
    echo "Computing for $TIME"
    date
    ${RUNCMD}$CODFILE  ${MODEL}_${TIME}_inter.in > ${MODEL}_${TIME}.out
    date
    #
    
    # Test if the month has finised properly
    echo "Test ${MODEL}_${TIME}.out"
    status=`tail -2 ${MODEL}_${TIME}.out | grep DONE | wc -l`
    if [[ $status == 1 ]]; then
      echo "All good"
      echo "XXXX${MYTEST}XXXX"
    else
      echo
      echo "Warning: month not finished properly"
      echo
      tail -20 ${MODEL}_${TIME}.out
      echo
      echo "Month ${TIME} did not work"
      echo
      exit 1
    fi
#
#  Archive
#
    LEVEL=0
    while [[ $LEVEL != $NLEVEL ]]; do
      if [[ ${LEVEL} == 0 ]]; then
        ENDF=
      else
        ENDF=.${LEVEL}
      fi
	  $CP -f ${MODEL}_rst.nc${ENDF} ${INIFILE}.nc${ENDF}
	  $MV -f ${MODEL}_his.nc${ENDF} ${MSSOUT}/${MODEL}_his_${TIME}.nc${ENDF}
	  $MV -f ${MODEL}_rst.nc${ENDF} ${MSSOUT}/${MODEL}_rst_${TIME}.nc${ENDF}
	  $MV -f ${MODEL}_avg.nc${ENDF} ${MSSOUT}/${MODEL}_avg_${TIME}.nc${ENDF}
      LEVEL=$((LEVEL + 1))
    done
    NM=$((NM + 1))
  done
  NY=$((NY + 1))
done
#
#############################################################












