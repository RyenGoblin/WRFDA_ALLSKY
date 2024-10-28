#!/bin/bash

export WRF_EM_CORE=1
echo 'WRF_EM_CORE=' $WRF_EM_CORE
export WRF_NMM_CORE=0
echo 'WRF_NMM_CORE=' $WRF_NMM_CORE
export WRF_DA_CORE=1
echo 'WRF_DA_CORE=' $WRF_DA_CORE
export OMP_NUM_THREADS=2
echo 'OMP_NUM_THREADS=' $OMP_NUM_THREADS
export OMP_STACKSIZE=2G
echo 'OMP_STACKSIZE=' $OMP_STACKSIZE
export PNETCDF_QUILT=1
echo 'PNETCDF_QUILT=' $PNETCDF_QUILT
ulimit -s unlimited
echo 'ulimit -s unlimited'



