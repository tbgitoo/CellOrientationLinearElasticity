from ngsolve import *
from ngsolve.bvp import BVP
import netgen.meshing as nm
import math
from numpy import linalg
from numpy import real
from os import chdir


# The idea here is that in the current directory,
# python can find the solution files u.sol, strain.sol etc.
def load_simulation(meshfile):
    m = nm.Mesh(3)
    m.Load(meshfile)
    m.SetBCName(0, "mid_chip_symmetry_plane_yz") # Face 1, bc index 1
    m.SetBCName(1, "gel_general_bottom_xy") # Face 2, bc index 2
    m.SetBCName(2, "end_gel_x_yz")  # Face 3, bc index 3
    m.SetBCName(3, "mid_channel_symmetry_plane_y_xz") # Face 4, bc index 4
    m.SetBCName(4, "gel_free_top_z_xy") # Face 5, bc index 5
    m.SetBCName(5, "mid_ridge_symmetry_plane_origin_xz") # Face 6, bc index 6 d
    m.SetBCName(6, "ridge_far_end_x_yz") # Face 7, bc index 7
    m.SetBCName(7, "rigde_long_vertical_face_y_xz") # Face 8, bc index 8
    m.SetBCName(8, "ridge_top_z_xy") # Face 9, bc index 9
    mesh = Mesh(m)
    fes = H1(mesh, order=4, dim=3)
    u=GridFunction(fes)
    u.Load("u.sol")
    fesp = H1(mesh, order=3, dim=9)
    strain = GridFunction(fesp)
    strain.Load("strain.sol")
    stretch = GridFunction(fesp)
    stretch.Load("stretch.sol")
    stress = GridFunction(fesp)
    stress.Load("stress.sol")
    return_dict = {
        "fes": fes,
        "u": u,
        "mesh": mesh,
        "strain": strain,
        "stretch": stretch,
        "stress": stress
    }
    return return_dict






