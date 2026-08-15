#-------------------------------------------------------------------------------
# XIOS
#-------------------------------------------------------------------------------

export USE_XIOS_ATM=0
export USE_XIOS_OCE=0

# Where to find the XIOS executable
export XIOS_EXE_DIR=${XIOS}/bin

# Process xml online
#-------------------
export ONLINE_XML="FALSE" 

# Names of xml files (defined in file_def_*)
#-------------------------------------------
export ATM_XIOS_NAME="wrfout wrf3d_1D wrf3d_1H" 
export OCE_OUTPUT_PREFIX="croco"
export OCE_XIOS_NAME="${OCE_OUTPUT_PREFIX}_3h_inst ${OCE_OUTPUT_PREFIX}_1h_avg_3d ${OCE_OUTPUT_PREFIX}_1h_inst_surf ${OCE_OUTPUT_PREFIX}_5d_aver" 

