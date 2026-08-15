############################# from mynamelist.sh ###############################

################################################################################
############################ USER CHANGES ######################################
################################################################################

# Run type (o/a/w, w.Afrc, oa, 2o1a, owa, owa.full...)
#  - Will select the models to use reading letters o/w/a/toywav/toyoce/toyatm 
#  - Will select the executables, and some options (see in the following sections)
#  - In coupled mode corresponds to the suffix of the OASIS_IN/namcouple.base.$RUNtype to use
export RUNtype=owa
export idealrun="FALSE" # set to true if it is an ideal run
export MOD=`echo $RUNtype | cut -d . -f 1`

#  Name of the experiment you are about to launch (max 30. CHAR)
export CEXPER=BENGUELA_${RUNtype}

#-------------------------------------------------------------------------------
# RUN_DIR
#-------------------------------------------------------------------------------

export EXEDIR_ROOT="$CWORK/rundir/${CEXPER}_execute"
export OUTPUTDIR_ROOT="$CWORK/rundir/${CEXPER}_outputs"
export RESTDIR_ROOT="$CWORK/rundir/${CEXPER}_restarts"

export  JOBDIR_ROOT=${CHOME}/jobs_${CEXPER}

