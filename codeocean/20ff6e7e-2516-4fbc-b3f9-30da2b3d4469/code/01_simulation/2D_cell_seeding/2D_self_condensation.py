from ngsolve import *
from ngsolve.bvp import BVP
import netgen.meshing as nm
import math
from os import chdir

import sys, os
sys.path.append('/code/01_simulation/2D_cell_seeding')
from simulation_cells_self_contraction import *

sys.path.append('/code/01_simulation')
from simulation_parameters import *



print("Simulation self-condensation 2D (effect of cell contraction)\n\n")
print("Parameters used\n")
print("Young modulus of the PDMS E=",E_PDMS,"\n")
print("Poission ratio of the PDMS nu=",nu_PDMS,"\n")
print("Young modulus of the cells E=",E_cells,"\n")
print("Poission ratio of the cells nu=",nu_cells,"\n")
print("Gravitational acceleration g_acc=",g_acc,"\n")
print("Isotropic cell contraction potential isotropic_cell_contraction=",isotropic_cell_contraction,"\n")
print("These simulations presume symmetry with walls remaing orthogonal to y at both y=0 and y=0.35. \n",
      "At y=0.35, y-displacement is possible, at y=0, we fix u=0\n")


cell_contraction=get_simulation_optimize_long_edge(rho=rho, E_PDMS=E_PDMS, nu_PDMS=nu_PDMS,
                                       E_cells=E_cells, nu_cells=nu_cells, g_acc=0,isotropic_cell_contraction=isotropic_cell_contraction)




cell_contraction["u"].Save("u.sol")
cell_contraction["stretch"].Save("stretch.sol")
cell_contraction["strain"].Save("strain.sol")
cell_contraction["stress"].Save("stress.sol")

#Draw(cell_contraction["u"],cell_contraction["mesh"],"u_self_contraction")



