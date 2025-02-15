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

n=4 # Power exponent for stretch distance addition

mesh=self_condensation["mesh"]
self_condensation_strain = GridFunction(H1(mesh, order=1, dim=9)) # keep it simple: order 1 means we can directly read
# the values for the gridpoints, no higher-order interpolation
self_condensation_strain.Set(self_condensation["strain"])
cyclic_chip_compression_0_strain = GridFunction(H1(mesh, order=1, dim=9)) # keep it simple: order 1 means we can directly read
# the values for the gridpoints, no higher-order interpolation
cyclic_chip_compression_0_strain.Set(cyclic_chip_compression_0["strain"])
cyclic_chip_compression_90_strain = GridFunction(H1(mesh, order=1, dim=9)) # keep it simple: order 1 means we can directly read
# the values for the gridpoints, no higher-order interpolation
cyclic_chip_compression_90_strain.Set(cyclic_chip_compression_90["strain"])


# Evauation of the in-plane (xy) direction that has the least and largest total compression



cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_0 = GridFunction(H1(mesh, order=1, dim=1))

most_compressive_angle_xy_self_condensation_and_cyclic_chip_compression_0 = GridFunction(H1(mesh, order=1, dim=1))

for ind in range(len(cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_0.vec)):
    cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_0.vec.data[ind]=least_compressive_direction_xy(
        cyclic_chip_compression_0_strain.vec[ind],self_condensation_strain.vec[ind],n=n
    )*180/math.pi
    most_compressive_angle_xy_self_condensation_and_cyclic_chip_compression_0.vec.data[ind]=most_compressive_direction_xy(
        cyclic_chip_compression_0_strain.vec[ind],self_condensation_strain.vec[ind],n=n
    )*180/math.pi

cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_90 = GridFunction(H1(mesh, order=1, dim=1))

most_compressive_angle_xy_self_condensation_and_cyclic_chip_compression_90 = GridFunction(H1(mesh, order=1, dim=1))

for ind in range(len(cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_90.vec)):
    cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_90.vec.data[
            ind] = least_compressive_direction_xy(
            cyclic_chip_compression_90_strain.vec[ind], self_condensation_strain.vec[ind], n=n
            ) * 180 / math.pi
    most_compressive_angle_xy_self_condensation_and_cyclic_chip_compression_90.vec[ind]=most_compressive_direction_xy(
        cyclic_chip_compression_90_strain.vec[ind],self_condensation_strain.vec[ind],n=n
    )*180/math.pi


vector_smallest_compression_xy_chip_compression_0=GridFunction(H1(mesh, order=1, dim=3))
for ind in range(len(vector_smallest_compression_xy_chip_compression_0.vec)):
    theAngle=cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_0.vec.data[ind]
    vector_smallest_compression_xy_chip_compression_0.vec.data[
            ind][0] = math.cos(-theAngle/180*math.pi)
    vector_smallest_compression_xy_chip_compression_0.vec.data[
        ind][1] = math.sin(-theAngle / 180 * math.pi)
    vector_smallest_compression_xy_chip_compression_0.vec.data[
        ind][2] = 0

vector_smallest_compression_xy_chip_compression_90 = GridFunction(H1(mesh, order=1, dim=3))
for ind in range(len(vector_smallest_compression_xy_chip_compression_90.vec)):
    theAngle = cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_90.vec.data[ind]
    vector_smallest_compression_xy_chip_compression_90.vec.data[
        ind][0] = math.cos(-theAngle / 180 * math.pi)
    vector_smallest_compression_xy_chip_compression_90.vec.data[
        ind][1] = math.sin(-theAngle / 180 * math.pi)
    vector_smallest_compression_xy_chip_compression_90.vec.data[
        ind][2] = 0



# stretch ratio
stretch_ratio_0= GridFunction(H1(mesh, order=1, dim=1))
for ind in range(len(stretch_ratio_0.vec)):
    strain_1=cyclic_chip_compression_0_strain.vec[ind]
    strain_2=self_condensation_strain.vec[ind]
    angle_min_compression=cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_0.vec[ind]/180*math.pi
    angle_max_compression=most_compressive_angle_xy_self_condensation_and_cyclic_chip_compression_0.vec[ind]/180*math.pi
    
    stretch_ratio_0.vec.data[ind]=stretch_ratio_xy(strain_1, strain_2,angle_min_compression,angle_max_compression,n=n)

# stretch ratio
stretch_ratio_90= GridFunction(H1(mesh, order=1, dim=1))
for ind in range(len(stretch_ratio_90.vec)):
    strain_1=cyclic_chip_compression_90_strain.vec[ind]
    strain_2=self_condensation_strain.vec[ind]
    angle_min_compression=cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_90.vec[ind]/180*math.pi
    angle_max_compression=most_compressive_angle_xy_self_condensation_and_cyclic_chip_compression_90.vec[ind]/180*math.pi
    
    stretch_ratio_90.vec.data[ind]=stretch_ratio_xy(strain_1, strain_2,angle_min_compression,angle_max_compression,n=n)


# Output mean orientation angles

unity_gridfunction=GridFunction(H1(mesh, order=1, dim=1))
unity_gridfunction.Set(1)
print("Mean orientation angle , 2D, cells on channel bottom, self-condensation and cyclic stretch in y-direction (perpendicular to grooves). ",
      "\n   Angle in degrees relative to stretch direction = ",
      90-Integrate(cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_90,mesh,definedon=mesh.Boundaries('cells_bottom_free_xy'))/
      Integrate(unity_gridfunction,mesh,definedon=mesh.Boundaries('cells_bottom_free_xy')),"\n\n")

print("Mean orientation angle, 2D, cells on channel bottom, self-condensation and cyclic stretch in x-direction (along grooves).",
      "\n   Angle in degrees relative to stretch direction = ",
      Integrate(cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_0,mesh,definedon=mesh.Boundaries('cells_bottom_free_xy'))/
      Integrate(unity_gridfunction,mesh,definedon=mesh.Boundaries('cells_bottom_free_xy')),"\n\n")







# Output to vtk

chdir("/results/evaluation/vtk")
vtk = VTKOutput(ma=mesh, coefs=[cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_0,most_compressive_angle_xy_self_condensation_and_cyclic_chip_compression_0, vector_smallest_compression_xy_chip_compression_0,stretch_ratio_0],
                names=["cell_orientation_angle","most_compressive_angle","cell_orientation_vector","stretch_ratio"], filename="2D_combined_condensation_stretch_0_degrees")

vtk.Do()

chdir("/results/evaluation/vtk")
vtk = VTKOutput(ma=mesh, coefs=[cell_orientation_angle_xy_self_condensation_and_cyclic_chip_compression_90, 
    most_compressive_angle_xy_self_condensation_and_cyclic_chip_compression_90,
    vector_smallest_compression_xy_chip_compression_90,stretch_ratio_90],
    names=["cell_orientation_angle","most_compressive_angle","cell_orientation_vector","stretch_ratio"], filename="2D_combined_condensation_stretch_90_degrees")

vtk.Do()


