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
from simulation_loading_3D import * 

chdir("/data/simulation_output/3D/cyclic_chip_compression_90_degrees")
cyclic_chip_compression_90=load_simulation("stretching_chip_mesh.vol")
chdir("/data/simulation_output/3D/self_condensation")
self_condensation=load_simulation("stretching_chip_mesh.vol")



# Evaluation for the combined effect of self-condensation and stretch
# ===================================================================

n=4 # Power exponent for stretch distance addition

mesh=self_condensation["mesh"]
self_condensation_strain = GridFunction(H1(mesh, order=1, dim=9)) # keep it simple: order 1 means we can directly read
# the values for the gridpoints, no higher-order interpolation
self_condensation_strain.Set(self_condensation["strain"])
cyclic_chip_compression_90_strain = GridFunction(H1(mesh, order=1, dim=9)) # keep it simple: order 1 means we can directly read
# the values for the gridpoints, no higher-order interpolation
cyclic_chip_compression_90_strain.Set(cyclic_chip_compression_90["strain"])



# Evauation of the in-plane (xy) direction that has the least and largest total compression



cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_90 = GridFunction(H1(mesh, order=1, dim=1))

most_compressive_angle_xy_self_condensation_and_cyclic_chip_compression_90 = GridFunction(H1(mesh, order=1, dim=1))

for ind in range(len(cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_90.vec)):
    cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_90.vec.data[ind]=least_compressive_direction_xy(
        cyclic_chip_compression_90_strain.vec[ind],self_condensation_strain.vec[ind],n=n
    )*180/math.pi
    most_compressive_angle_xy_self_condensation_and_cyclic_chip_compression_90.vec.data[ind]=most_compressive_direction_xy(
        cyclic_chip_compression_90_strain.vec[ind],self_condensation_strain.vec[ind],n=n
    )*180/math.pi



# From the angle, also produce the orientation vector for the arrows
vector_smallest_compression_xy_chip_compression_90=GridFunction(H1(mesh, order=1, dim=3))
for ind in range(len(vector_smallest_compression_xy_chip_compression_90.vec)):
    theAngle=cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_90.vec.data[ind]
    vector_smallest_compression_xy_chip_compression_90.vec.data[
            ind][0] = math.cos(-theAngle/180*math.pi)
    vector_smallest_compression_xy_chip_compression_90.vec.data[
        ind][1] = math.sin(-theAngle / 180 * math.pi)
    vector_smallest_compression_xy_chip_compression_90.vec.data[
        ind][2] = 0

# stretch ratio
stretch_ratio_90= GridFunction(H1(mesh, order=1, dim=1))
for ind in range(len(stretch_ratio_90.vec)):
    strain_1=cyclic_chip_compression_90_strain.vec[ind]
    strain_2=self_condensation_strain.vec[ind]
    angle_min_compression=cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_90.vec[ind]/180*math.pi
    angle_max_compression=most_compressive_angle_xy_self_condensation_and_cyclic_chip_compression_90.vec[ind]/180*math.pi
    
    stretch_ratio_90.vec.data[ind]=stretch_ratio_xy(strain_1, strain_2,angle_min_compression,angle_max_compression,n=n)

# Output mean orientation angle

channel_gel_indicator=IfPos(y-0.175,1,0)*IfPos(5-x,1,0)*IfPos(0.1-z,1,0)

unity_gridfunction=GridFunction(H1(mesh, order=1, dim=1))
unity_gridfunction.Set(1)
print("Mean orientation angle , 3D, cells in gel in grooves, lower 100 micron layer, self-condensation and cyclic stretch in y-direction (perpendicular to the grooves). ",
      "\n   Angle in degrees relative to stretch direction = ",
      90-Integrate(cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_90*channel_gel_indicator,mesh)/
      Integrate(channel_gel_indicator,mesh),"\n\n")


# Output to vtk

chdir("/results/evaluation/vtk")
vtk = VTKOutput(ma=mesh, coefs=[cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_90,most_compressive_angle_xy_self_condensation_and_cyclic_chip_compression_90, vector_smallest_compression_xy_chip_compression_90,stretch_ratio_90],
                names=["cell_orientation_angle","most_compressive_angle","cell_orientation_vector","stretch_ratio"], filename="3D_combined_condensation_stretch_90_degrees")

vtk.Do()




