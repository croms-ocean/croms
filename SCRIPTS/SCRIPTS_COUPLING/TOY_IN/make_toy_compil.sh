#!/bin/bash

source ../myenv_mypath.sh
ln -sf Makefile.${MACHINE} Makefile

for toytype in atm oce wav ; do
    ln -sf read_namelist.F90_${toytype} read_namelist.F90
    ln -sf toy_model.F90_${toytype} toy_model.F90
    make clean > toy_${toytype}_clean.out
    make > toy_${toytype}_make.out
    mv toy_model toy_${toytype}
done

