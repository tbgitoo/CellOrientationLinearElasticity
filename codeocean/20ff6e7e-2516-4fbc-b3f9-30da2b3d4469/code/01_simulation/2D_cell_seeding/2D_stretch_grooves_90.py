from ngsolve import *
from ngsolve.bvp import BVP
import netgen.meshing as nm
import math
from os import chdir

import sys, os
sys.path.append('/code/01_simulation/2D_cell_seeding')
from simulation_cells import *


sys.path.append('/code/01_simulation')
from simulation_parameters import *


print("Simulation externally applied stretch to cells seed '2D' on to the PDMS\n\n")
print("Parameters used\n")
print("Young modulus of the PDMS E=",E_PDMS,"\n")
print("Poisson ratio of the PDMS nu=",nu_PDMS,"\n")
print("Young modulus of the cells E=",E_cells,"\n")
print("Poisson ratio of the cells nu=",nu_cells,"\n")
print("Gravitational acceleration g_acc=",g_acc,"\n")
print("Chip compression ratio chip_compression=",chip_compression,"\n")
print("Orientation: Chip compression perpendicular to the groove direction\n")

cyclic_chip_compression=get_simulation(rho=1000, E_PDMS=E_PDMS, nu_PDMS=nu_PDMS,
                                       E_cells=E_cells, nu_cells=nu_cells, g_acc=g_acc,chip_compression=chip_compression,isotropic_cell_contraction=0,
                   ridge_lateral_partial_adhesion=True,orientation_stretch=math.pi/4)


cyclic_chip_compression["u"].Save("u.sol")
cyclic_chip_compression["stretch"].Save("stretch.sol")
cyclic_chip_compression["strain"].Save("strain.sol")
cyclic_chip_compression["stress"].Save("stress.sol")




