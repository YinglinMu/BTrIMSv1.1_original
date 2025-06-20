#!/bin/bash
#PBS -N wvfac_PBL
#PBS -P w40
#PBS -q normalsr
#PBS -l walltime=20:00:00
#PBS -l ncpus=104
#PBS -l mem=350GB
#PBS -l jobfs=100MB
#PBS -l storage=gdata/hh5+gdata/w28+gdata/rt52+gdata/w40+gdata/w97+scratch/w40
#PBS -l wd

ulimit -s unlimited 
module load intel-compiler netcdf
export OMP_STACKSIZE=2G
export OMP_NUM_THREADS=104

rm -f *.o *.mod
# compile all files
ifort -c global_data.f90
ifort -c util.f90 
ifort -c input_data_handling_era5.f90 
ifort -c bt_subs.f90 
ifort -c -qopenmp back_traj.f90 

ifort -qopenmp *.o -lnetcdff -lnetcdf -o main 

#./main 07 10 2023 08 10 2023 /scratch/w40/ym7079/test0_PBL/ > terminal_output.txt


