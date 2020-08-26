"""
testing stub to import dockeyeMS_energy module and 
use energy_c method
"""
import dockeyeMS_energy
n1 = 5
n2 = 6
nmod = 7
atom_data = [1.,2.,3.]
#atom_data = 5.
x = dockeyeMS_energy.energy_c(n1,n2,nmod,atom_data)
print(type(x))
print(x)
