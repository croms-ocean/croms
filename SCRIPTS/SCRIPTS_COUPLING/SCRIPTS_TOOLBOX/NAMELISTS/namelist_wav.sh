#-------------------------------------------------------------------------------
# WAV
#-------------------------------------------------------------------------------

# Where to find the wav executable
if [[ $RUNtype =~ .*owa.* ]] ; then
    export WAV_EXE_DIR=${WAV}/exe_owa
elif [[ $RUNtype =~ .*ow.* ]] ; then
    export WAV_EXE_DIR=${WAV}/exe_ow
elif [[ $RUNtype =~ .*aw.* ]] ; then
    export WAV_EXE_DIR=${WAV}/exe_aw
else
    export WAV_EXE_DIR=${WAV}/exe_frc
fi

# Namelist
#---------
# Chosing the ww3_shel.inp.base.SHELL_EXT (see options in WW3_IN)
if [[ $RUNtype =~ .*toy.* ]] ; then
    export SHELL_EXT=$MOD
else
    export SHELL_EXT=$RUNtype
fi

# Time steps
#-----------
export DT_WAV=3600     # TMAX = 3*TCFL
export DT_WW_PRO=1200  # TCFL = 0.8 x dx/(g/fmin4pi) with fmin=0.0373 => 3-4 % of dx
export DT_WW_REF=1800  # TMAX / 2
export DT_WW_SRC=10    # TSRC = usually 10s  (could be between 5s and 60s)

# Grid size
#----------
export wavnx=41 ; export wavny=42

# Parameter
#---------- 
export hmin=75; # e.g. minimum water depth in CROCO (will be used to delimit coastline in WW3)

# Forcing files
#--------------
# forcin: forcing file(s) PREFIX list (input file are supposed to be in the form: PREFIX_Y????M??.nc)
# forcww3: name of ww3_prnc.inp extension, e.g current or wind/era5, see in WW3_IN directory
if [[ $RUNtype =~ .*owa.* ]]; then
    export forcin=()
    export forcww3=()
elif [[ $RUNtype =~ .*Afrc.* || $RUNtype =~ .*ow.* ]] ; then
    export forcin=(ERA5_wind)
    export forcww3=(wind.era5)
elif [[ $RUNtype =~ .*Ofrc.* ]] ; then
    export forcin=(CROCO_current CROCO_level)
    export forcww3=(current level)
elif [[ $RUNtype =~ .*frc.* ]] ; then
    export forcin=(ERA5_wind CROCO_current CROCO_level)
    export forcww3=(wind.era5 current level)
fi

# Boundary files
#---------------
# prefix for boundary files (leave empty is none), there are supposed to be in WAV_FILES_DIR
export bouncin= 

# Output settings
#----------------
export wav_int=21600  # output interval (s)
export wav_pnt=0  # point output interval. Put 0 if no point output
export point_output_list=${WAV_FILES_DIR}/my_point_output_test.txt # file where to find list of point (format: lon lat name) to output spectrum
export wav_trck=0     # track output interval. Put 0 if no track output  
export flagout="TRUE" # Keep (TRUE) or not (FALSE) the ww3 output binary files (e.g. out_grd.ww3)

