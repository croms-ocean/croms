#!/bin/bash
#
################### USER'S DEFINED PARAMETERS ##########################
#
#=======================================================================
#  Define environment variables
#=======================================================================
#
# Environment for Crontab
#
OS=`uname`
echo "OPERATING SYSTEM IS: $OS"

if [ -f /etc/bashrc ]; then
  . /etc/bashrc
fi
if [[ $OS == Linux ]] ; then
  source ~/.bashrc
elif [[ $OS == Darwin ]] ; then
  source ~/.bash_profile
else
  echo "unknown operating system" 1>&2
fi
#
export TOOLSDIR=$HOME/croco/croco_tools/Forecast_tools
export RUNDIR=${PWD} 
export MATLAB=/usr/local/bin/matlab
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/loaddap-3.5.2/lib

unalias cp
unalias mv
export CP=/bin/cp
export MV=/bin/mv
export LN=/bin/ln
export DATESTR=`date +%Y-%m-%d`
#
#=======================================================================
#  Define options
#=======================================================================
#
echo Start Forecast 
date
#
# Clean results of previous forecasts
#
export CLEAN=0
#
# Get forcing files from DODS SERVER and process them for CROCO
# PRE_PROCESS=1 ==> do the work (0 otherwise)
#
export PRE_PROCESS=1
#
# Restart from previous forecast
#
export RESTART=0
#
# Run hindcast/forecast
#
export RUN=1
#
#
# Make a few plots
export PLOT=1
#
#=======================================================================
# Define time parameters (it will modify croco_forecast.in)
#=======================================================================
#
# Model time step [seconds]
DT=3600
# Number of barotropic time steps within one baroclinic time step [number]
# NDTFAST in croco.in
NFAST=60
# Hindcast depth [days] see crocotools_param (hdays/fdays)
NDAYS_HIND=1
# Forecast depth [days]
NDAYS_FCST=3
# Output frequency [hours]
#   average
ND_AVG=24
#   history (if = -1 set equal to NUMTIMES)
ND_HIS=6
#   restart (if = -1 set equal to NUMTIMES)
ND_RST=24
# 
#=======================================================================
#  Define directories, files and running parameters
#=======================================================================
#
export SCRATCHDIR=${RUNDIR}/SCRATCH
export INPUTDIR=${RUNDIR}
export MSSDIR=${RUNDIR}/CROCO_FILES
export MSSOUT=${RUNDIR}/FORECAST
#
export MODEL=croco
export CODFILE=./croco
export EXEC="mpirun -np 4 "
#
export GRDFILE=${MODEL}_grd.nc
export INIFILE=${MODEL}_ini.nc
export RSTFILE=${MODEL}_rst.nc
export AVGFILE=${MODEL}_avg.nc
export HISFILE=${MODEL}_his.nc
export BLKFILE=${MODEL}_blk_GFS_0.nc
export FRCFILE=${MODEL}_frc_GFS_0.nc
export BRYFILE=${MODEL}_bry_mercator_0.nc
export CLMFILE=${MODEL}_clm_mercator_0.nc
export INIFILE0=${MODEL}_ini_mercator_0.nc
#
################  END USER'S DEFINED PARAMETERS ##########################
#
# Go the Input directory 
# 
cd $INPUTDIR
#
# Create working directories if needed
#
if [[ ! -e $SCRATCHDIR ]]; then
    mkdir $SCRATCHDIR
else
    echo "$SCRATCHDIR already exists" 1>&2
fi
if [[ ! -e $MSSOUT ]]; then
    mkdir $MSSOUT
else
    echo "$MSSOUT already exists" 1>&2
fi
#
# Clean I/O files from previous forecast
#
if [ $CLEAN = 1 ] ; then
  rm -f $MSSDIR/${MODEL}_blk_* $MSSDIR/${MODEL}_frc_* $MSSDIR/${MODEL}_bry_*
  rm -f $MSSDIR/${MODEL}_clm_* $MSSDIR/${MODEL}_ini_*
  rm -f $MSSOUT/${MODEL}_his_* $MSSOUT/${MODEL}_avg_* $MSSOUT/${MODEL}_rst_*
fi
#
#=======================================================================
# PRODUCE INITIAL CONDITIONS AND FORCING
#=======================================================================
# Compute lateral boundaries from MERCATOR
#     and surface forcing from GFS
#
if [ $PRE_PROCESS = 1 ]; then
  echo "Processing boundary and forcing files"
  [ ! -e start.m ] && { echo "==> Stop : you need a start.m file"; exit 1 ; }
  $MATLAB -nodisplay -batch "run('start.m') ; run('make_forecast.m');" > matlab_forecast.out
fi
#
# Copy files in SCRATCH dir
#
echo "Getting $BLKFILE from $MSSDIR"
$CP -f  $MSSDIR/$BLKFILE $SCRATCHDIR
echo "Getting $FRCFILE from $MSSDIR"
$CP -f  $MSSDIR/$FRCFILE $SCRATCHDIR
echo "Getting $BRYFILE from $MSSDIR"
$CP -f  $MSSDIR/$BRYFILE $SCRATCHDIR
echo "Getting $CLMFILE from $MSSDIR"
$CP -f  $MSSDIR/$CLMFILE $SCRATCHDIR
#
# Initial file:
#  if RESTART, get initial conditions 
#  from the previous hindcast run
#  else use global data.
#
if [ $RESTART = 1 ] ; then
  echo "Getting $INIFILE from $MSSOUT"
  $CP -f  $MSSOUT/$INIFILE $SCRATCHDIR
else
  echo "Getting $INIFILE from $MSSDIR"
  $CP -f  $MSSDIR/$INIFILE0 $SCRATCHDIR/$INIFILE
fi
#
# Get static files
#
echo "Getting $CODFILE from $INPUTDIR"
$CP -f $INPUTDIR/$CODFILE $SCRATCHDIR
chmod u+x $CODFILE
echo "Getting ${GRDFILE} from $MSSDIR"
$CP -f $MSSDIR/${GRDFILE} $SCRATCHDIR
echo "Getting ${MODEL}_forecast.in from $INPUTDIR"
$CP -f $INPUTDIR/${MODEL}_forecast.in $SCRATCHDIR
echo "Getting ${MODEL}_stations.in from $INPUTDIR"
$CP -f $INPUTDIR/${MODEL}_stations.in $SCRATCHDIR
#
# Time management in croco_forecast.in
#
NDAYS=$((NDAYS_HIND + NDAYS_FCST ))
NUMTIMES=$((NDAYS * 24 * 3600 / DT))
[[ ${ND_AVG} -ne -1 ]] &&  NUMAVG=$((ND_AVG * 3600 / DT )) || NUMAVG=$NUMTIMES
[[ ${ND_HIS} -ne -1 ]] &&  NUMHIS=$((ND_HIS * 3600 / DT )) || NUMHIS=$NUMTIMES
[[ ${ND_RST} -ne -1 ]] &&  NUMRST=$((ND_RST * 3600 / DT )) || NUMRST=$NUMTIMES

sed -e 's/NUMTIMES/'$NUMTIMES'/' -e 's/TIMESTEP/'$DT'/' -e 's/NFAST/'$NFAST'/' \
    -e 's/NUMAVG/'$NUMAVG'/' -e 's/NUMHIS/'$NUMHIS'/' -e 's/NUMRST/'$NUMRST'/' \
    -e 's/<logfilename>/'${MODEL}_${TIME}.out'/' \
        < $INPUTDIR/${MODEL}_forecast.in > $SCRATCHDIR/${MODEL}_forecast.in
#
#  Change directory
#
cd $SCRATCHDIR
#
#=======================================================================
#  HINDCAST/FORECAST RUN
#=======================================================================
#
if [ $RUN = 1 ] ;then
  echo Hindcast/Forecast run
  date
  $EXEC $CODFILE ${MODEL}_forecast.in > ${MODEL}_forecast_${DATESTR}.out
  date
  #
  # Store the initial file for the next forecast
  #
  $CP -f $SCRATCHDIR/$RSTFILE ${MSSOUT}/$INIFILE
  #
  # Store the output forecast files
  #
  $CP $SCRATCHDIR/$HISFILE ${MSSOUT}/${MODEL}_his_forecast_${DATESTR}.nc
  $CP $SCRATCHDIR/$AVGFILE ${MSSOUT}/${MODEL}_avg_forecast_${DATESTR}.nc
fi
#
#=======================================================================
#  PLOT RESULTS
#=======================================================================
#
if [ $PLOT = 1 ] ;then
  #
  cd $INPUTDIR
  [ ! -e start.m ] && { echo "==> Stop : you need a start.m file"; exit 1 ; }
  # 
  # Production plot
  [ ! -e plot_forecast_croco.m ] && $CP -f ${TOOLSDIR}/plot_forecast_croco.m .
  $MATLAB -nodisplay -batch  "run('plot_forecast_croco.m');" > plot_forecast_croco.out
  #
fi

echo Forecast finished
date

