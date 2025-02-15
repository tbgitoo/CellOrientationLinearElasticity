from ngsolve import *
from ngsolve.bvp import BVP
import netgen.meshing as nm
import math
from numpy import linalg
from numpy import real
from os import chdir
from CellOrientationLinearElasticity import *

import sys, os
sys.path.append('/code/02_evaluation')
from simulation_loading_2D import * 


chdir("/data/simulation_output/2D/cyclic_chip_compression_0_degrees")
cyclic_chip_compression_0=load_simulation("cells_2D_seeded.vol")
chdir("/data/simulation_output/2D/cyclic_chip_compression_90_degrees")
cyclic_chip_compression_90=load_simulation("cells_2D_seeded.vol")
chdir("/data/simulation_output/2D/self_condensation")
self_condensation=load_simulation("cells_2D_seeded.vol")



# Evaluation of orientation for cyclic chip compression only

mesh=self_condensation["mesh"]
s=self_condensation["strain"]
s_order_1 = GridFunction(H1(mesh, order=1, dim=9))
s_order_1.Set(s)
strain_ratio_self_condensation = GridFunction(H1(mesh, order=1, dim=1))
eigenvector_smallest_compression_xy_self_condensation=GridFunction(H1(mesh, order=1, dim=3))
eigenvector_smallest_compression_self_condensation=GridFunction(H1(mesh, order=1, dim=3))

for ind in range(len(strain_ratio_self_condensation.vec)):
    strain_ratio_self_condensation.vec.data[ind]=eigenvalue_ratio_xy(s_order_1.vec[ind])
    m = least_compressive_eigenvector_xy(s_order_1.vec[ind])
    if(abs(m[0]) > abs(m[1]) and real(m[0])<0):
        m[0]=-m[0]
        m[1]=-m[1]
    if (abs(m[1]) > abs(m[0]) and real(m[1]) < 0):
        m[0] = -m[0]
        m[1] = -m[1]
    for ind_dim in range(len(m)):
        eigenvector_smallest_compression_xy_self_condensation.vec.data[ind][ind_dim]=-real(m[ind_dim])
    m = least_compressive_eigenvector(s_order_1.vec[ind])
    if (abs(m[0]) > abs(m[1]) and real(m[0]) < 0):
        m[0] = -m[0]
        m[1] = -m[1]
        m[2] = -m[2]
    if (abs(m[1]) > abs(m[0]) and real(m[1]) < 0):
        m[0] = -m[0]
        m[1] = -m[1]
        m[2] = -m[2]
    for ind_dim in range(len(m)):
        eigenvector_smallest_compression_self_condensation.vec.data[ind][ind_dim] = -real(m[ind_dim])




# Cell orientation in degrees by atan from the xy vector for self-condensation
cell_orientation_angle_self = GridFunction(H1(mesh, order=1, dim=1))
for ind in range(len(cell_orientation_angle_self.vec)):
    mx=abs(eigenvector_smallest_compression_xy_self_condensation.vec[ind][0])
    my = abs(eigenvector_smallest_compression_xy_self_condensation.vec[ind][1])
    if mx==0:
        cell_orientation_angle_self.vec.data[ind]=90
    else:
        cell_orientation_angle_self.vec.data[ind] =  math.atan(my/mx)/math.pi*180




# Evaluate orientation in cyclic compression along the grooves (orientation 0)

mesh=cyclic_chip_compression_0["mesh"]
s=cyclic_chip_compression_0["strain"]
s_order_1 = GridFunction(H1(mesh, order=1, dim=9))
s_order_1.Set(s)
strain_ratio_chip_compression_0 = GridFunction(H1(mesh, order=1, dim=1))
eigenvector_smallest_compression_xy_chip_compression_0=GridFunction(H1(mesh, order=1, dim=3))
eigenvector_smallest_compression_chip_compression_0=GridFunction(H1(mesh, order=1, dim=3))

for ind in range(len(strain_ratio_chip_compression_0.vec)):
    strain_ratio_chip_compression_0.vec.data[ind]=eigenvalue_ratio_xy(s_order_1.vec[ind])
    m = least_compressive_eigenvector_xy(s_order_1.vec[ind])
    if (abs(m[0]) > abs(m[1]) and real(m[0]) < 0):
        m[0] = -m[0]
        m[1] = -m[1]
    if (abs(m[1]) > abs(m[0]) and real(m[1]) < 0):
        m[0] = -m[0]
        m[1] = -m[1]
    for ind_dim in range(len(m)):
        eigenvector_smallest_compression_xy_chip_compression_0.vec.data[ind][ind_dim]=-real(m[ind_dim])
    m = least_compressive_eigenvector(s_order_1.vec[ind])
    if (abs(m[0]) > abs(m[1]) and real(m[0]) < 0):
        m[0] = -m[0]
        m[1] = -m[1]
        m[2] = -m[2]
    if (abs(m[1]) > abs(m[0]) and real(m[1]) < 0):
        m[0] = -m[0]
        m[1] = -m[1]
        m[2] = -m[2]
    for ind_dim in range(len(m)):
        eigenvector_smallest_compression_chip_compression_0.vec.data[ind][ind_dim] = -real(m[ind_dim])



cell_orientation_angle_0 = GridFunction(H1(mesh, order=1, dim=1))
for ind in range(len(cell_orientation_angle_0.vec)):
    mx=abs(eigenvector_smallest_compression_xy_chip_compression_0.vec[ind][0])
    my = abs(eigenvector_smallest_compression_xy_chip_compression_0.vec[ind][1])
    if mx==0:
        cell_orientation_angle_0.vec.data[ind]=90
    else:
        cell_orientation_angle_0.vec.data[ind] =  math.atan(my/mx)/math.pi*180




# Evaluate orientation in cyclic compression perpendicular to the grooves (orientation 90)

mesh=cyclic_chip_compression_90["mesh"]
s=cyclic_chip_compression_90["strain"]
s_order_1 = GridFunction(H1(mesh, order=1, dim=9))
s_order_1.Set(s)
strain_ratio_chip_compression_90 = GridFunction(H1(mesh, order=1, dim=1))

# Evaluate also using the same algorithm as for
# the combined strain approach
strain_ratio_chip_compression_90_combination_algorithm = GridFunction(H1(mesh, order=1, dim=1))

eigenvector_smallest_compression_xy_chip_compression_90=GridFunction(H1(mesh, order=1, dim=3))
eigenvector_smallest_compression_chip_compression_90=GridFunction(H1(mesh, order=1, dim=3))

for ind in range(len(strain_ratio_chip_compression_90.vec)):
    strain_ratio_chip_compression_90.vec.data[ind]=eigenvalue_ratio_xy(s_order_1.vec[ind])
    m = least_compressive_eigenvector_xy(s_order_1.vec[ind])
    if (abs(m[0]) > abs(m[1]) and real(m[0]) < 0):
        m[0] = -m[0]
        m[1] = -m[1]
    if (abs(m[1]) > abs(m[0]) and real(m[1]) < 0):
        m[0] = -m[0]
        m[1] = -m[1]
    for ind_dim in range(len(m)):
        eigenvector_smallest_compression_xy_chip_compression_90.vec.data[ind][ind_dim]=-real(m[ind_dim])
    m = least_compressive_eigenvector(s_order_1.vec[ind])
    if (abs(m[0]) > abs(m[1]) and real(m[0]) < 0):
        m[0] = -m[0]
        m[1] = -m[1]
        m[2] = -m[2]
    if (abs(m[1]) > abs(m[0]) and real(m[1]) < 0):
        m[0] = -m[0]
        m[1] = -m[1]
        m[2] = -m[2]
    for ind_dim in range(len(m)):
        eigenvector_smallest_compression_chip_compression_90.vec.data[ind][ind_dim] = -real(m[ind_dim])



cell_orientation_angle_90 = GridFunction(H1(mesh, order=1, dim=1))
for ind in range(len(cell_orientation_angle_90.vec)):
    mx=abs(eigenvector_smallest_compression_xy_chip_compression_90.vec[ind][0])
    my = abs(eigenvector_smallest_compression_xy_chip_compression_90.vec[ind][1])
    if mx==0:
        cell_orientation_angle_90.vec.data[ind]=90
    else:
        cell_orientation_angle_90.vec.data[ind] =  math.atan(my/mx)/math.pi*180





# Output to vtk

chdir("/results/evaluation/vtk")
vtk = VTKOutput(ma=mesh, coefs=[cyclic_chip_compression_0["u"],cell_orientation_angle_0,
                eigenvector_smallest_compression_xy_chip_compression_0, eigenvector_smallest_compression_chip_compression_0,
                strain_ratio_chip_compression_0                ],
                names=["u","cell_orientation_angle","preferred_cell_orientation_xy","preferred_cell_orientation", "strain_ratio_xy"
                       ], filename="2D_cyclic_chip_compression_0")

vtk.Do()

chdir("/results/evaluation/vtk")
vtk = VTKOutput(ma=mesh, coefs=[cyclic_chip_compression_90["u"],cell_orientation_angle_90,
                eigenvector_smallest_compression_xy_chip_compression_90, eigenvector_smallest_compression_chip_compression_90,
                strain_ratio_chip_compression_90                ],
                names=["u","cell_orientation_angle","preferred_cell_orientation_xy","preferred_cell_orientation", "strain_ratio_xy"
                       ], filename="2D_cyclic_chip_compression_90")


vtk.Do()



chdir("/results/evaluation/vtk")
vtk = VTKOutput(ma=mesh, coefs=[self_condensation["u"],cell_orientation_angle_self,
                eigenvector_smallest_compression_xy_self_condensation, eigenvector_smallest_compression_self_condensation,
                strain_ratio_self_condensation                ],
                names=["u","cell_orientation_angle","preferred_cell_orientation_xy","preferred_cell_orientation", "strain_ratio_xy"
                       ], filename="2D_self_condensation")


vtk.Do()





