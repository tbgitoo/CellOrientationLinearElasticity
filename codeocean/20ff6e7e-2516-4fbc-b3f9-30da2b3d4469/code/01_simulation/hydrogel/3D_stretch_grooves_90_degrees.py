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

print("Simulation externally applied stretch 3D hydrogel geometry\n\n")
print("Parameters used\n")
print("Young modulus of the dECM-fibrin hydrogel E=",E_hydrogel,"\n")
print("Long-term Poission ratio of the dECM-fibrin hydrogel nu=",nu_hydrogel_self_condensation,"\n")
print("Gravitational acceleration g_acc=",g_acc,"\n")
print("Chip compression ratio chip_compression=",chip_compression,"\n")
print("Orientation: Chip compression perpendicular to the groove direction\n")


cyclic_chip_compression=get_simulation(rho=1000, E=E_hydrogel, nu=nu_hydrogel_chip_compression, g_acc=g_acc,
                                       chip_compression=chip_compression, isotropic_cell_contraction=0,
                                       ridge_lateral_partial_adhesion=False,orientation_stretch=math.pi/4)


cyclic_chip_compression["u"].Save("u.sol")
cyclic_chip_compression["stretch"].Save("stretch.sol")
cyclic_chip_compression["strain"].Save("strain.sol")
cyclic_chip_compression["stress"].Save("stress.sol")




