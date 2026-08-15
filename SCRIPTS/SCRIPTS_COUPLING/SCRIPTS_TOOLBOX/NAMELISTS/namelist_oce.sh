#-------------------------------------------------------------------------------
# OCE
#-------------------------------------------------------------------------------

# Where to find or put the croco execuatble
export OCE_EXE_DIR=${CHOME}/CROCO_IN

# Online Compilation
#-------------------
#!!!!!!! IMPORTANT NOTE !!!!!!!
# If activated, creates croco executable depending on this namelist. 
#   - In param.h it modifies the grid size, the number of procs in x and y direction with those given in myjob.sh
#   - In cppdefs.h it modifies the following options with informations given below
#       MPI, OA_COUPLING, OW_COUPLING, MRL_WCI, 
#       XIOS, LOGFILE, MPI_NOLAND,
#       AGRIF, AGRIF_2WAY, 
#       BULK_FLUX, ONLINE, AROME, ARPEGE, ERA_ECMWF
#       FRC_BRY, CLIMATOLOGY
#       TIDES, PSOURCE, PSOURCE_NCFILE, PSOURCE_NCFILE_TS
# Other changes of parameterizations, numerical schemes, etc should be made "by hand" in CROCO_IN/cppdefs.h.base
# NB: ONLINE_COMP=1 does not work when RESTART_FLAG="TRUE" in myjob.sh
#     (this case may arise when running a coupled run after a forced spinup for instance)
#     in such case, the user should compile "by hand" prior to the run
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
export ONLINE_COMP=1

# Time steps
#-----------
export DT_OCE=3600
export NDTFAST=60

# EXACT_RESTART
#--------------
export EXACT_RESTART="FALSE"

# MPI_NOLAND
#-----------
export MPI_NOLAND="FALSE"

# Domains
#--------
if [[ $RUNtype =~ .*2o.* ]] ; then
# Simulation with 2 ocean domains
    export AGRIFZ=1
    export AGRIF_2WAY="TRUE"
else
# Simulation with 1 ocean domain
    export AGRIFZ=0
    export AGRIF_2WAY="FALSE"
fi

# Atm coupling
#--------------
if [[ $RUNtype =~ .*NOpsfc.* ]] ; then
    export CPL_PATM="FALSE"
else
    export CPL_PATM="TRUE"
fi

# Wave coupling
#--------------
# OW_COUPLING_FULL
#   - TRUE: uses all terms from the wave model (stokes, tke flux, bernouilli head, orbital velocities)
#   - FALSE: uses their monochromatic approximation
if [[ $RUNtype =~ .*full.* ]] ; then
    export OW_COUPLING_FULL="TRUE"
else
    export OW_COUPLING_FULL="FALSE"
fi
# WAVE_SMFLUX
#   - TRUE: uses stress from the wave model 
#   - FALSE: uses stress from the bulk formulation
#   Not used when cpl with atm as the stress comes from the atm model. For owa coupled cases both atm and wav stresses are used. 
if [[ $RUNtype =~ .*two.* ]] ; then
    export WAVE_SMFLUX="TRUE"
else
    export WAVE_SMFLUX="FALSE"
fi

# Forcings
#---------
export ini_ext='ini_SODA' # ini extension file (ini_SODA,...)
export bdy_ext='bry_SODA' # bry extension file (clm_SODA,bry_SODA,...)

# flag for surface forcing should be true except in the case of atm coupling
if [[ $MOD =~ .*a.* ]] ; then
    export surfrc_flag="FALSE"
else
    export surfrc_flag="TRUE"
fi
export interponline=0 # switch (1=on, 0=off) for online surface interpolation. Only works with MONTHLY input files!
export frc_ext='blk_ERA5' # surface forcing extension(blk_ERA5, frc_ERA5,...). If interponline=1 precise the type (ERA_ECMWF or AROME,  [CFSR by default], names as cppkey name in croco)

export tide_flag="FALSE" # the forcing extension must be blk_??? otherwise tide forcing overwrites it 
export river_flag="FALSE"

# Output settings
#----------------
#!!! WARNING: when XIOS is activated the following values (for the model) are not taken into account
export oce_his_sec=86400     # history output interval (in number of second) 
export oce_avg_sec=86400     # average output interval (in number of second) 

