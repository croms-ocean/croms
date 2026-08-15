#-------------------------------------------------------------------------------
# TOY
#-------------------------------------------------------------------------------

# Where to find the toy executable(s)
export TOY_EXE_DIR=${CHOME}/TOY_IN

# Choose for which model you use the toy 
# If several separate with spaces
# options are: oce atm wav
#---------------------------------------
export toytype=("wav" "atm") 

# Forcing files that will be read by the toy
#-------------------------------------------
export toyfile=("$CWORK/TOY_FILES/ww3.200501.nc" \
                "$CWORK/TOY_FILES/wrfout_d01_2005-01-01_00:00:00")
export toytimerange=('2,10' \
                     '2,10')

