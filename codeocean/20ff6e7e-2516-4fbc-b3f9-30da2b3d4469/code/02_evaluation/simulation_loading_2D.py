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
    m.SetBCName(0, "origin_PDMS_yz") # Face 1
    m.SetBCName(1, "origin_PDMS_xz") # Face 2
    m.SetBCName(2, "interior_1_between_PDMS_and_cell_bottom_yz")
    m.SetBCName(3, "ridge_PDMS_yz")
    m.SetBCName(4, "cells_bottom_free_vertical_xz")
    m.SetBCName(5, "inner_2_between_PDMS_and_cell_bottomx_xy")
    m.SetBCName(6, "PDMS_bottom_free_surface_xz")
    m.SetBCName(7, "cells_bottom_free_xy")
    m.SetBCName(8, "cells_bottom_symmetry_yz")
    m.SetBCName(9, "cells_bottom_mid_channel_xz")
    m.SetBCName(10, "cells_top_far_end_yz")
    m.SetBCName(11, "cells_top_y_end_xz")
    m.SetBCName(12, "cells_top_free_xy")
    m.SetBCName(13, "cells_top_far_end_to_ridge_xz")
    m.SetBCName(14, "cells_top_ridge_end_yz")
    m.SetBCName(15, "cells_top_ridge_end_xz")
    m.SetBCName(16, "cells_top_symmetry_yz")
    m.SetBCName(17, "PDMS_far_end_yz")
    m.SetBCName(18, "PDMS_y_end_xz")
    m.SetBCName(19, "inner_3_between_PDMS_and_cell_top_xy")
    m.SetBCName(20, "PDMS_bottom_xy")
    m.SetBCName(21, "ridge_PDMS_xz")
    m.SetMaterial(1,"PDMS")
    m.SetMaterial(2, "cells_top")
    m.SetMaterial(3, "cells_bottom")
    mesh = Mesh(m)
    fes = H1(mesh, order=2, dim=3,dgjumps=True)
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






