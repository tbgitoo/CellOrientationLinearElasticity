from ngsolve import *
from ngsolve.bvp import BVP
import netgen.meshing as nm
import math
from os import chdir


import sys, os
sys.path.append('/code/01_simulation/hydrogel')
from simulation_hydrogel import *

sys.path.append('/code/01_simulation')
from simulation_parameters import *

print("Simulation self-condensation 3D\n\n")
print("Parameters used\n")
print("Young modulus of the dECM-fibrin hydrogel E=",E_hydrogel,"\n")
print("Long-term Poission ratio of the dECM-fibrin hydrogel nu=",nu_hydrogel_self_condensation,"\n")
print("Gravitational acceleration g_acc=",g_acc,"\n")
print("Isotropic cell contraction potential isotropic_cell_contraction=",isotropic_cell_contraction,"\n")



cell_contraction=get_simulation(rho=rho, E=E_hydrogel, nu=nu_hydrogel_self_condensation, g_acc=g_acc,
                                chip_compression=0, isotropic_cell_contraction=isotropic_cell_contraction,
                                       ridge_lateral_partial_adhesion=True,orientation_stretch=0)

cell_contraction["u"].Save("u.sol")
cell_contraction["stretch"].Save("stretch.sol")
cell_contraction["strain"].Save("strain.sol")
cell_contraction["stress"].Save("stress.sol")

#Draw(cell_contraction["u"],cell_contraction["mesh"],"u_self_contraction")



