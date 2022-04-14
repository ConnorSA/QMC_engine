## @package energy_plotting

## @brief This file just outputs the energy against the bondlength

## Imported modules
import numpy as np
import matplotlib.pyplot as plt

## Loading in the data from the energies.txt
## @param data contains the bond lengths and energy conatined in the energies.txt file
data = np.loadtxt('energies.txt')
print('energies.txt uploaded successfully')

## @param bond_length contains all the bond lengths evaluated 
bond_length = data[:,0]
## @param energy conatins the energy results evaluated at their corresponding bond_lengths
energy = data[:,1]

## @brief creates the plot for bond length on the x axis and the resulting energy on the y axis
plt.figure()
plt.plot(bond_length, energy)
plt.xlabel('Bond Length (Atomic Units)'); plt.ylabel('Energy (Hartree units)')
plt.title('Energy of 2 atom system with varying bond lengths')
plt.show()
