#!/bin/bash 

module swap PrgEnv-cray PrgEnv-gnu
module load cray-netcdf-hdf5parallel/4.4.1
module load cray-hdf5-parallel/1.10.0
module load fftw
module load cray-petsc-complex 
module swap cray-mpich/7.0.4 cray-mpich/7.4.4


