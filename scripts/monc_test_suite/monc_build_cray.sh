#!/bin/bash 

compiler=$1
system=$2

echo $compiler

if [ $system = "MO-XC40" ]; then

    if [ $compiler = "cray" ]; then
	#source set_env_cray_cray.sh
	module load cray-netcdf-hdf5parallel/4.4.1
	module load cray-hdf5-parallel/1.10.0
	module load fftw
	module load cray-petsc-complex
	module swap cce/8.3.4 cce/8.4.3
	module swap cray-mpich/7.0.4 cray-mpich/7.4.4
	module list
	fcm make -j32 -f fcm-make/monc-cray-cray.cfg -f fcm-make/casim_socrates_mirror.cfg --new
    else 
	#source set_env_cray_gnu.sh
#	source /common/lfric/modules/setup
#	module load meto-environment/lfric/gnu/6.2.0
#	module load fftw
#	module load cray-petsc-complex
	module swap PrgEnv-cray PrgEnv-gnu
	module load cray-netcdf-hdf5parallel/4.4.1
	module load cray-hdf5-parallel/1.10.0
	module load fftw
	module load cray-petsc-complex 
	module swap cray-mpich/7.0.4 cray-mpich/7.4.4
	module list
	fcm make -j32 -f fcm-make/monc-cray-gnu.cfg -f fcm-make/casim_socrates_mirror.cfg --new
    fi
fi

