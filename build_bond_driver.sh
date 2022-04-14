#!/bin/bash
#chmod +rx build_script.sh #if you have permission issues.
#List here all the f90 files to compile, separated by spaces
#THERE MUST NOT be spaces around '=' in bash
#Driver goes last in the list.

#Starts python gui to get user input and generate txt file
#Not sure if this goes here or within the makefile
# python3 init_params.py

myprogramfiles="constants.f90 component_functions.f90 basis_functions.f90 param_search.f90 MCMC.f90  read_data.f90 electron_density.f90 calc.f90 netcdf_file.f90 bond_driver.f90"
# 
#Name of compiled file
outfile="bond_driver_exec"

#Name of compiler
fc=gfortran
#Use nf-config to grab the compile and link flags. Backticks run command and grab output
fflags=`nf-config --fflags`
flibs=`nf-config --flibs`
#Actual compile line. Other flags etc can be added
$fc  -g $fflags $myprogramfiles $flibs -fopenmp -O2 -o $outfile  -std=f2008 # -Wall

#mv *.mod clean_up/.
