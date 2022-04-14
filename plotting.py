## @package plotting

## @brief Main plotting script which outputs the contour plots for the wavefunction or electron probability density

## Modules to be loaded
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from netCDF4 import Dataset

## Uploading the xy plane coordinates from the xyz.txt file generated from the xy_grid subroutine in netcdf_file.f90
coords = np.loadtxt('xyz.txt')

## Uploading the main results NetCDF file named results.nc4 generated from the result_netcdf routine
## @param data this contains all the information contained in the NetCDF file
data = nc.Dataset('results.nc4', mode='r', format='NETCDF4')
## Prints a statement if it has been uploaded successfully
print('File result.nc4 uploaded successfully')
print(' ')
print(data)
print(' ')

## @brief Extracting the number of electrons and nuclei variables from the result file

## @param num_ele is the number of electrons used
num_ele = data.variables['Num_of_Electrons'][:]
print('Number of electrons: ', num_ele)
## @param num_nuc is the number of nuclei
num_nuc = data.variables['Num_of_Nuclei'][:]
print('Number of nuclei: ', num_nuc)


## Extracts the x and y coordinates from the xyz.txt
## @param x is the x coordinates
## @param y is the y coordinates
x = coords[:,0]
y = coords[:,1]

## @brief This does the contour plot for the 1 electron case. 
## It returns a contour plot of the wavefunction squared in the xy plane
if num_ele == 1.0:
    
    ## @brief Extracts the wavefunction and degrees of freedom parameters from the input NetCDF file
    ## @param wavefunction contains the wavefunction data
    wavefunction = data.variables['Electron_Density'][:]
    ## @param dof contains the degrees of freedom data
    dof = data.variables['Optimal_DOF'][:]
    print('Optimal Degrees of freedom', dof)

    ## Creates the contour plot of the wavefunction squared
    plt.figure()
    plt.tricontourf(x, y, wavefunction**2)
    plt.colorbar()
    plt.xlabel('X (Atomic Units)'); plt.ylabel('Y (Atomic Units)')
    plt.title('Electron Probability Density ')
    plt.show()
    
## @brief This does the contour plot for the 2 electron case.  
## It returns a contour plot of the electron probability density in the xy plane   
elif num_ele == 2.0:
    
    ## @param ele_den the electron probability density
    ele_den = data.variables['Electron_Density'][:]
    ## @param dof the degrees of freedom
    dof = data.variables['Optimal_DOF'][:]

    ## Creates the contour plot of the electron probability density
    plt.figure()
    plt.tricontourf(x, y, ele_den)
    plt.colorbar()
    plt.xlabel('X (Atomic Units)'); plt.ylabel('Y (Atomic Units)')
    plt.title('Electron Probability Density ')
    plt.show()
    
## If there is an issue with the number of electrons given, this error statement will print
else:
    print('Number of electrons not defined')
    print('Error in plotting')


