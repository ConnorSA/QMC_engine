#!/bin/bash
#chmod +rx build_script.sh #if you have permission issues.
#List here all the f90 files to compile, separated by spaces
#THERE MUST NOT be spaces around '=' in bash
#Driver goes last in the list.

#Starts python gui to get user input and generate txt file
#Not sure if this goes here or within the makefile
# python3 init_params.py

myprogramfiles="constants.f90 component_functions.f90 basis_functions.f90 param_search.f90 MCMC.f90 read_data.f90 GP_surrogate.f90 gradient_estimator.f90 stoch_grad.f90  Biased_Optim.f90 Biased_Optim_example_driver.f90  "
#"latin_driver.f90"

#Name of compiled file
outfile="Bias_optim_driver_exec"

#Name of compiler
fc=gfortran
#Use nf-config to grab the compile and link flags. Backticks run command and grab output, O2 allows for compiler optimisation
fflags=`nf-config --fflags`
# Includes openmp and mkl libraries, will need to change path to where mkl is located  
flibs=`nf-config --flibs` 
#Actual compile line. Other flags etc can be added
$fc  -g $fflags $myprogramfiles $flibs -fopenmp -o $outfile -std=f2008 -llapack

#mv *.mod clean_up/.