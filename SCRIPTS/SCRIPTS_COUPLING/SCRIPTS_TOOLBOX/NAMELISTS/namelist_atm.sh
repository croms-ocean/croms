#-------------------------------------------------------------------------------
# ATM
#-------------------------------------------------------------------------------

# Where to find the atm exectuable
if [[ $RUNtype =~ .*a.* && ( $RUNtype =~ .*o.* || $RUNtype =~ .*w.* ) ]] ; then
    export ATM_EXE_DIR=${ATM}/exe_coupled
else
    export ATM_EXE_DIR=${ATM}/exe_uncoupled
fi

# Namelist
#---------
#!!!!!!! IMPORTANT NOTE !!!!!!!
# Changes of parameterizations, numerical schemes, etc in atmnamelist
# should be made "by hand" in the WRF_IN/namelist.input.base.complete file
# Only settings into <...> in WRF_IN/namelist.input.base.complete are automatically
# filled in by the present mynamelist settings
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
export atmnamelist=namelist.input.base.complete

# Time steps
#-----------
export DT_ATM=150

# Boundaries 
#-----------
export interval_seconds=21600 # interval ( in sec ) of the lateral input
export auxinput4_interval=360 # interval ( in min ) of bottom input
export nbmetsoil=4
export nbmetlevel=38

# Physics
#--------
 # Cd formulation (default = 0, wave cpl needs = 5)
if [[ $RUNtype =~ .*aw.* || $RUNtype =~ .*owa.* ]] ; then
    export isftcflx=5
else
    export isftcflx=1
fi

# Domains
#--------
export NB_dom=1 # Number of coupled domains
export wrfcpldom='d01' # which WRF domain to couple
export nestfeedback="TRUE" # 1 way (FALSE) or 2 Way (TRUE) nesting
export onlinecplmask="TRUE" # Erase existing CPLMASK and build default mask (depending on the nb of atm and oce domains)

# Moving nest 
#------------
export ATM_CASE="DEFAULT"  # no moving nest: DEFAULT or with: MOVING_NEST
# if ATM_CASE=DEFAULT, the following is not used
export num_mv_nest=0 # number of moving nests
# if several nest, the following variables need to have the format "1st_nest 2nd_nest"
export ref_coef="3" # refinement coef for nest
export ew_size="283" # nest size in east-west dim ([multiple of ref_coef] + 1)
export ns_size="295" # nest size in north-south dim ([multiple of ref_coef] + 1)
export i_prt_strt="580" # where nest is starting in parent's grid x-dim 
export j_prt_strt="59" # where nest is starting in parent's grid y-dim
# Tracking parameters
export vortex_interval=5 # When to update vortex position
export max_vortex_speed=40 # Used to compute the search radius for the new vortex center position
export corral_dist=8 # The closest distance between child and parend boundary (in parent grid cell)
export track_level=50000 # The pressure level (in Pa) where the vortex is tracked
export time_to_move=0 # The time (in minutes) until the nest is moved (at the beginning)

# Nudging (assimilation) options
#-------------------------------
export switch_fdda=0 # To activate fdda nudging
export nudgedom="1" # select which kind of nudging you want (1=grid-nudging, 2=spectral nudging) for each domain. Example for spectral nudging over parent only "2 0"
export nudge_coef="0.0003" # nudge coef. Need to be the same size than nudge
export nudge_interval_m="360" # time interval (in min) between analysis times 
export nudge_end_h="144" # time (in hours) to stop nudging after start of forecast

# Output settings
#----------------
#!!! WARNING: when XIOS is activated the following values (for the model) are not taken into account
export atm_his="180 180 60"               # output interval (minutes) for each domain
export atm_his_frames=1000 # $((31*24))   # nb of outputs per file
export atm_diag_int_m=${atm_his}          # diag output interval (m)
export atm_diag_frames=1000               # nb of diag outputs per file
# file for specifying different than default output variables: OPTIONAL, leave empty if not used
export atm_iofields='myoutfields.txt'

