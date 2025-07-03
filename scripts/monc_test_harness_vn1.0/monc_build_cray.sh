#!/bin/bash 

compiler=$1
system=$2
opt_level=$3
modules_load=$4
compiler_dir=$5

echo ''
echo "Building MONC with:"
echo "  system       : $system"
echo "  compiler     : $compiler"
echo "  modules      : $modules_load"
echo "  optimisation : $opt_level"
echo ''

### load the modules for appropriate compiler

if [ $system = "MO-XC40" ]; then

    if [ $compiler = "cray" ]; then
        module swap PrgEnv-gnu PrgEnv-cray
	if [ $modules_load = "default" ]; then 
	    module cray-netcdf-hdf5parallel/4.4.1
	    module load cray-hdf5-parallel/1.10.0
	    module load fftw
	    module load cray-petsc-complex
	    module swap cce/8.3.4 cce/8.4.3
	    module swap cray-mpich/7.0.4 cray-mpich/7.4.4
            module load perftools
	else
	    module load cdt/18.12
	    module load cray-netcdf-hdf5parallel
	    module load cray-hdf5-parallel
	    module load fftw
	    module load cray-petsc-complex
	    module load cray-tpsl 
            module load perftools-base
            module load perftools
	fi
        module list
        echo ""
	### Build depending on debug level
	if [ $opt_level = 'debug' ]; then
	    echo ' Building cray in debug (no optimisation). This will be slow.'
	    fcm make -j36 -f fcm-make/monc-cray-cray-debug.cfg -f fcm-make/casim_socrates.cfg --new
	elif [ $opt_level = 'safe' ]; then 
	    echo ' Building cray in safe. This will be slow but should be vector safe (optimisation = 2).'
	    fcm make -j36 -f fcm-make/monc-cray-cray-safe.cfg -f fcm-make/casim_socrates.cfg --new
	elif [ $opt_level = 'high' ]; then
	    echo ' Building cray in high. This will be fast but not vector safe (optimisation = 3).'
	    fcm make -j36 -f fcm-make/monc-cray-cray.cfg -f fcm-make/casim_socrates.cfg --new
	else
	    echo ' No debug level has been set. Model will not build and harness will not work.'
	fi
        echo ""

    elif [ $compiler = "intel" ]; then
	### no module_loads condition as cdt intel modules do not load properly
	### and intel does not work with monc
	module swap PrgEnv-cray PrgEnv-intel
	module swap intel/15.0.0.090 intel/18.0.5.274
	module load cray-netcdf-hdf5parallel/4.4.1
	module load cray-hdf5-parallel/1.10.0
	module load fftw
	module load cray-petsc-complex 
	module swap cray-mpich/7.0.4 cray-mpich/7.4.4
        module list
        echo ""
	if [ $opt_level = 'debug' ]; then
            echo ' Building intel in debug (no optimisation). This will be slow.'
	    fcm make -j36 -f fcm-make/monc-cray-intel-debug.cfg -f fcm-make/casim_socrates.cfg --new
	elif [ $opt_level = 'safe' ]; then 
            echo ' Building intel in safe. This will be slow but should be vector safe (optimisation = 2).'
	    fcm make -j36 -f fcm-make/monc-cray-intel-safe.cfg -f fcm-make/casim_socrates.cfg --new
	elif [ $opt_level = 'high' ]; then
            echo ' Building intel in high. This will be fast but not vector safe (optimisation = 3).'
	    fcm make -j36 -f fcm-make/monc-cray-intel.cfg -f fcm-make/casim_socrates.cfg --new
	else
            echo ' No debug level has been set. Model will not build and harness will not work.'
	fi
        echo ""

    else #gnu 
        if [ $modules_load = "default" ]; then 
	    module swap PrgEnv-cray PrgEnv-gnu
	    module load cray-netcdf-hdf5parallel/4.4.1
	    module load cray-hdf5-parallel/1.10.0
	    module load fftw
	    module load cray-petsc-complex 
	    module swap cray-mpich/7.0.4 cray-mpich/7.4.4
	else
	    module load cdt/18.12
	    module swap PrgEnv-cray PrgEnv-gnu
	    module load cray-netcdf-hdf5parallel
	    module load cray-hdf5-parallel
	    module load fftw
	    module load cray-petsc-complex
	    module load cray-tpsl  
	fi    
        module list
        echo ""
	### Build depending on debug level
	if [ $opt_level = 'debug' ]; then
            echo ' Building gnu in debug (no optimisation). This will be slow.'
	    fcm make -j36 -f fcm-make/monc-cray-gnu-debug.cfg -f fcm-make/casim_socrates.cfg --new
	elif [ $opt_level = 'safe' ]; then 
            echo ' Building gnu in safe. This will be slow but should be vector safe (optimisation = 2).'
	    fcm make -j36 -f fcm-make/monc-cray-gnu-safe.cfg -f fcm-make/casim_socrates.cfg --new
	elif [ $opt_level = 'high' ]; then
            echo ' Building gnu in high. This will be fast but not vector safe (optimisation = 3).' 
	    fcm make -j36 -f fcm-make/monc-cray-gnu.cfg -f fcm-make/casim_socrates.cfg --new
	else
            echo ' No debug level has been set. Model will not build and harness will not work.'
	fi
        echo ""

    fi

elif [ $system = "ARCHER2" ]; then

    if [ $compiler = "cray" ]; then
        if [ $modules_load = "default" ]; then
            source env/archer_cray_mod
        else
            echo " only default modules are configured for ARCHER2."
            exit 1
        fi # modules check
        module list
        echo ""
        ### Build depending on debug level
        if [ $opt_level = 'debug' ]; then
            echo ' Building cray in debug (no optimisation). This will be slow.'
            fcm make -j4 -f fcm-make/monc-cray-cray-debug.cfg -f fcm-make/casim_socrates.cfg --new
        elif [ $opt_level = 'safe' ]; then
            echo ' Building cray in safe. This will be slow but should be vector safe (optimisation = 2).'
            fcm make -j4 -f fcm-make/monc-cray-cray-safe.cfg -f fcm-make/casim_socrates.cfg --new
        elif [ $opt_level = 'high' ]; then
            echo ' Building cray in high. This will be fast but not vector safe (optimisation = 3).'
            fcm make -j4 -f fcm-make/monc-cray-cray.cfg -f fcm-make/casim_socrates.cfg --new
        else
            echo ' No debug level has been set. Model will not build and harness will not work.'
        fi
        echo ""

    else # gnu
        if [ $modules_load = "default" ]; then
            source env/archer_gnu_mod
        else
            echo " only default modules are configured for ARCHER2."
            exit 1
        fi # modules check
        module list
        echo ""

        ### Build depending on debug level
        if [ $opt_level = 'debug' ]; then
            echo ' Building gnu in debug (no optimisation). This will be slow.'
            fcm make -j4 -f fcm-make/monc-cray-gnu-debug.cfg -f fcm-make/casim_socrates.cfg --new
        elif [ $opt_level = 'safe' ]; then
            echo ' Building gnu in safe. This will be slow but should be vector safe (optimisation = 2).'
            fcm make -j4 -f fcm-make/monc-cray-gnu-safe.cfg -f fcm-make/casim_socrates.cfg --new
        elif [ $opt_level = 'high' ]; then
            echo ' Building gnu in high. This will be fast but not vector safe (optimisation = 3).' 
            fcm make -j4 -f fcm-make/monc-cray-gnu.cfg -f fcm-make/casim_socrates.cfg --new
        else
            echo ' No debug level has been set. Model will not build and harness will not work.'
        fi # debug check
        echo ""

    fi # compiler check
fi # system check

echo ''

# Build the compilation record:
date > version
echo "" >> version
echo "Built on the '${system}' system with the '${compiler}' compiler, optimisation level: '${opt_level}'" >> version
echo "" >> version
echo "Loaded modules with '${modules_load}' option:"
(module list) 2>> version
echo "" >> version
echo "fcm info:" >> version
fcm info >> version
echo "" >> version
echo "fcm status:" >> version
fcm status >> version
echo "" >> version
echo "fcm diff:" >> version
fcm diff >> version


