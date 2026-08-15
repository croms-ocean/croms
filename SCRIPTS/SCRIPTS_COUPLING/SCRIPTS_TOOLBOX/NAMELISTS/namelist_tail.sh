################################################################################
############################ END USER CHANGE ###################################
################################################################################

#-------------------------------------------------------------------------------
# Setting models used
#-------------------------------------------------------------------------------
export MOD=`echo $RUNtype | cut -d . -f 1`
if [[ $RUNtype =~ .*toywav.* ]] ; then
    export USE_TOYWAV=1
    export USE_WAV=0
elif [[ $MOD =~ .*w.* ]] ; then
    export USE_TOYWAV=0
    export USE_WAV=1
else
    export USE_TOYWAV=0
    export USE_WAV=0
fi

if [[ $RUNtype =~ .*toyatm.* ]] ; then
    export USE_TOYATM=1
    export USE_ATM=0
elif [[ $MOD =~ .*a.* ]] ; then
    export USE_TOYATM=0
    export USE_ATM=1
else
    export USE_TOYATM=0
    export USE_ATM=0
fi

if [[ $RUNtype =~ .*toyoce.* ]] ; then
    export USE_TOYOCE=1
    export USE_OCE=0
elif [[ $MOD =~ .*o.* ]] ; then
    export USE_TOYOCE=0
    export USE_OCE=1
else
    export USE_TOYOCE=0
    export USE_OCE=0
fi

# Avoid error when a model is not used
[[ -z ${USE_ATM+x} ]] && export USE_ATM=0
[[ -z ${USE_OCE+x} ]] && export USE_OCE=0
[[ -z ${USE_WAV+x} ]] && export USE_WAV=0
[[ -z ${USE_XIOS_ATM+x} ]] && export USE_XIOS_ATM=0
[[ -z ${USE_XIOS_OCE+x} ]] && export USE_XIOS_OCE=0
[[ -z ${USE_TOYATM+x} ]] && export USE_TOYATM=0
[[ -z ${USE_TOYOCE+x} ]] && export USE_TOYOCE=0
[[ -z ${USE_TOYWAV+x} ]] && export USE_TOYWAV=0
[[ -z ${ONLINE_XML+x} ]] && export ONLINE_XML="FALSE"

# KEY for XIOS and TOY #
export USE_XIOS=$(( ${USE_XIOS_ATM} + ${USE_XIOS_OCE} ))
export USE_TOY=$(( ${USE_TOYATM} + ${USE_TOYOCE} + ${USE_TOYWAV} ))
[ ${USE_TOY} -ge 1 ] && export USE_CPL=1 || export USE_CPL=$(( ${USE_ATM} * ${USE_OCE} + ${USE_ATM} * ${USE_WAV} + ${USE_OCE} * ${USE_WAV} ))

echo '  '
echo 'Model(s) used:'
for use in USE_ATM USE_OCE USE_WAV USE_XIOS_ATM USE_XIOS_OCE USE_TOYATM USE_TOYOCE USE_TOYWAV ; do 
    if [[ $use -ge 1 ]] ; then 
        echo $use 
    fi
done

if [ ${USE_ATM} == 0 ]; then
    export DT_ATM=1
fi

if [ ${USE_OCE} == 0 ]; then
    export DT_OCE=1
fi

if [ ${USE_WAV} == 0 ]; then
    export DT_WAV=1
fi

#-------------------------------------------------------------------------------
# Setting toy stuff
#-------------------------------------------------------------------------------

[ -z ${toytype+x} ] && export toytype=()
export nbtoy=${#toytype[@]}
export model_to_toy=()
export toynamelist=()
export DT_TOY=()

if [[ ${USE_TOY} -ge 1 ]]; then
    [[ ${USE_TOY} != $nbtoy ]] && { echo "There is an incoherence between the number of USE_TOY_* activated and the number of arguement in toytype. Make Sure they are coherent." ; exit 1;}
    for k in `seq 0 $(( ${nbtoy} - 1))` ; do
        [ ${toytype[$k]} == "oce" ] && model_to_toy+=("croco")
        [ ${toytype[$k]} == "atm" ] && model_to_toy+=("wrf")
        [ ${toytype[$k]} == "wav" ] && model_to_toy+=("ww3")
        toynamelist+=("TOYNAMELIST.nam.${toytype[$k]}.${MOD}")
        if [ ${toytype[$k]} == "oce" ]; then
            targ=$( ncdump -v time ${toyfile[$k]} | grep "time = " | sed -n '2p' | cut -d ' ' -f 4-5 )
            tsp=$(( $( echo "${targ}" | cut -d ',' -f 2) - $( echo "${targ}" | cut -d ',' -f 1) ))
            DT_TOY+=("${tsp}")
        elif [ ${toytype[$k]} == "atm" ]; then
            targ=$( ncks -v XTIME ${toyfile[$k]} | grep "XTIME ="| cut -d '=' -f 2 | cut -d ',' -f 1-2)
            tsp=$(( $( echo "${targ}" | cut -d ',' -f 2)*60 - $( echo "${targ}" | cut -d ',' -f 1)*60 ))
            DT_TOY+=("${tsp}")
        elif [ ${toytype[$k]} == "wav" ]; then
            targ=$( ncdump -v time ${toyfile[$k]} | grep "time =" | sed -n '2p' | cut -d ' ' -f 4-5)
            arg1=$( echo "${targ}" | cut -d ',' -f 2 )
            arg2=$( echo "${targ}" | cut -d ',' -f 1 )
            tsp=$( echo "(${arg1}  - ${arg2})*86400" | bc | cut -d '.' -f 1)
            DT_TOY+=("${tsp}")
        fi
    done
fi

#-------------------------------------------------------------------------------
# Checking exp. options
#-------------------------------------------------------------------------------
echo "  "
echo " Checking experiment options... "

if [ ${USE_CPL} -ge 1 ]; then
    if [ $(( ${CPL_FREQ} % ${DT_ATM} )) -ne 0 ] || \
       [ $(( ${CPL_FREQ} % ${DT_OCE} )) -ne 0 ] || \
       [ $(( ${CPL_FREQ} % ${DT_WAV} )) -ne 0 ] ; then
        printf "\n\n ERROR: Problem of consistency between Coupling Frequency and Time Step with ATM, OCE or WAV model, we stop...\n\n" && exit 1
    fi
    if [ ${USE_TOY} -ge 1 ]; then
        for k in `seq 0 $(( ${nbtoy} - 1))` ; do
            if [ $(( ${CPL_FREQ} % ${DT_TOY[$k]} )) -ne 0 ] ; then
                printf "\n\n ERROR: Problem of consistency between Coupling Frequency and Time Step for TOY model, we stop...\n\n" && exit 1
            fi
        done
    fi
fi

