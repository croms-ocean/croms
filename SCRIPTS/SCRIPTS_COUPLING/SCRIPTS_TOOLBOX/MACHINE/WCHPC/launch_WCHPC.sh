############ CREATE app.conf for DATARMOR ###########

if [ ${USE_ATM} -eq 1 ]; then
    echo "-np ${NP_ATM} ./wrfexe" >> app.conf
    if [ ${USE_XIOS_ATM} -eq 1 ]; then
        echo "-np ${NP_XIOS_ATM} ./xios_server.exe" >> app.conf
    fi
fi

if [ ${USE_OCE} -eq 1 ]; then
    if [[ ${MPI_NOLAND} == "TRUE" ]]; then
        echo "-np ${NP_OCE} ./crocox croco.in" >> app.conf
    else
        echo "-np $(( ${NP_OCEX} * ${NP_OCEY} )) ./crocox croco.in" >> app.conf
    fi

    if [ ${USE_XIOS_OCE} -eq 1 ]; then
        echo "-np ${NP_XIOS_OCE} ./xios_server.exe" >> app.conf
    fi
fi

if [ ${USE_WAV} -eq 1 ]; then
    echo "-np ${NP_WAV} ./wwatch" >> app.conf
fi

if [ ${USE_TOY} -ge 1 ]; then
    for k in `seq 0 $(( ${nbtoy} - 1 ))`; do
        echo "-np ${NP_TOY} ./toy${toytype[$k]}" >> app.conf     
    done
fi
##
