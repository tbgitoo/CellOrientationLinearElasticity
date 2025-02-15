from ngsolve import *
from ngsolve.bvp import BVP
import netgen.meshing as nm
import math
from CellOrientationLinearElasticity import *



def get_simulation_with_y_on_vertical_long_face(rho=1000, E_PDMS=2e6, nu_PDMS=0.49, E_cells=3e4, nu_cells=0.49, g_acc=9.81,
                                                chip_compression=0.2,isotropic_cell_contraction=0.05,
                                                ridge_lateral_partial_adhesion=True, orientation_stretch=0,y_on_long_edge=0):
    m = nm.Mesh(3)
    m.Load("/code/01_simulation/2D_cell_seeding/cells_2D_seeded.vol")
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
    domain_values_E = {'PDMS': E_PDMS, 'cells_top': E_cells, 'cells_bottom': E_cells}
    E = CoefficientFunction([domain_values_E[mat]
                   for mat in mesh.GetMaterials()])
    domain_values_nu = {'PDMS': nu_PDMS, 'cells_top': nu_cells, 'cells_bottom': nu_cells}
    nu = CoefficientFunction([domain_values_nu[mat]
                             for mat in mesh.GetMaterials()])
    mu = E / 2 / (1 + nu)
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
     # We use here the penalty dictionnary schemes to fix the relevant coordinates
    fes = H1(mesh, order=2, dim=3,dgjumps=True)
    g = Parameter(g_acc)
    density = Parameter(rho)
    # Gravity pull
    coef_force_z = CoefficientFunction((0, 0, -density * g))
    v = fes.TestFunction()
    u = fes.TrialFunction()
    f = LinearForm(fes)
    f += coef_force_z * v * dx
    a=getElasticBilinearForm(fes=fes, E=E, nu=nu)
    # Penalties: only hold within planes on the bottom and symmetry. Due to cell contraction, the lateral surfaces
    # of the ridges can detach (y penalty of inner_origin_xz and inner_y_xz
    penalty_x_dict = {"origin_PDMS_yz": 1e20,
                      "origin_PDMS_xz": 0,
                      "interior_1_between_PDMS_and_cell_bottom_yz": 0,
                      "ridge_PDMS_yz": 0,
                      "cells_bottom_free_vertical_xz": 0,
                      "inner_2_between_PDMS_and_cell_bottomx_xy": 0,
                      "PDMS_bottom_free_surface_xz": 0,
                      "cells_bottom_free_xy": 0,
                      "cells_bottom_symmetry_yz": 1e20,
                      "cells_bottom_mid_channel_xz": 0,
                      "cells_top_far_end_yz":1e20,
                      "cells_top_y_end_xz":0,
                      "cells_top_free_xy":0,
                      "cells_top_far_end_to_ridge_xz":0,
                      "cells_top_ridge_end_yz":0,
                      "cells_top_ridge_end_xz":0,
                      "cells_top_symmetry_yz":1e20,
                      "PDMS_far_end_yz":1e20,
                      "PDMS_y_end_xz":0,
                      "inner_3_between_PDMS_and_cell_top_xy":0,
                      "PDMS_bottom_xy":0,
                      "ridge_PDMS_xz": 0
                      }
    penalty_y_dict = {"origin_PDMS_yz": 0,
                      "origin_PDMS_xz": 1e20,
                      "interior_1_between_PDMS_and_cell_bottom_yz": 0,
                      "ridge_PDMS_yz": 0,
                      "cells_bottom_free_vertical_xz": 0,
                      "inner_2_between_PDMS_and_cell_bottomx_xy": 0,
                      "PDMS_bottom_free_surface_xz": 0,
                      "cells_bottom_free_xy": 0,
                      "cells_bottom_symmetry_yz": 0,
                      "cells_bottom_mid_channel_xz": 1e20,
                      "cells_top_far_end_yz": 0,
                      "cells_top_y_end_xz": 0,
                      "cells_top_free_xy": 0,
                      "cells_top_far_end_to_ridge_xz": 1e20,
                      "cells_top_ridge_end_yz": 0,
                      "cells_top_ridge_end_xz": 0,
                      "cells_top_symmetry_yz": 0,
                      "PDMS_far_end_yz": 0,
                      "PDMS_y_end_xz": 0,
                      "inner_3_between_PDMS_and_cell_top_xy": 0,
                      "PDMS_bottom_xy": 0,
                      "ridge_PDMS_xz": 0
                      }
    penalty_z_dict = {"origin_PDMS_yz": 0,
                      "origin_PDMS_xz": 0,
                      "interior_1_between_PDMS_and_cell_bottom_yz": 0,
                      "ridge_PDMS_yz": 0,
                      "cells_bottom_free_vertical_xz": 0,
                      "inner_2_between_PDMS_and_cell_bottomx_xy": 0,
                      "PDMS_bottom_free_surface_xz": 0,
                      "cells_bottom_free_xy": 0,
                      "cells_bottom_symmetry_yz": 0,
                      "cells_bottom_mid_channel_xz": 0,
                      "cells_top_far_end_yz": 0,
                      "cells_top_y_end_xz": 0,
                      "cells_top_free_xy": 0,
                      "cells_top_far_end_to_ridge_xz": 0,
                      "cells_top_ridge_end_yz": 0,
                      "cells_top_ridge_end_xz": 0,
                      "cells_top_symmetry_yz": 0,
                      "PDMS_far_end_yz": 0,
                      "PDMS_y_end_xz": 0,
                      "inner_3_between_PDMS_and_cell_top_xy": 0,
                      "PDMS_bottom_xy": 1e20,
                      "ridge_PDMS_xz": 0
                      }
    penalty_x = CoefficientFunction([penalty_x_dict[bound] for bound in mesh.GetBoundaries()])
    penalty_y = CoefficientFunction([penalty_y_dict[bound] for bound in mesh.GetBoundaries()])
    penalty_z = CoefficientFunction([penalty_z_dict[bound] for bound in mesh.GetBoundaries()])
    a += penalty_x * u[0] * v[0] * ds
    a += penalty_y * u[1] * v[1] * ds
    a += penalty_z * u[2] * v[2] * ds
    # Cyclic chip compression: we indicate an angle to the main groove direction (x) here. We here apply the
    # boundary value addition for target displacement as a function of the penalties valid at the given boundary
    if orientation_stretch==math.pi/4:
        f += penalty_y * (-chip_compression)*y* v[1] * ds
    else:
        if orientation_stretch==0:
            f += penalty_x * (-chip_compression) * x * v[0] * ds
        else:
            f += penalty_x * (-chip_compression*math.cos(orientation_stretch)) * x * v[0] * ds
            f += penalty_y * (-chip_compression * math.sin(orientation_stretch)) * y * v[1] * ds
    domain_values_E_cells_only = {'PDMS': 0, 'cells_top': E_cells, 'cells_bottom': E_cells}
    E_cells_only = CoefficientFunction([domain_values_E_cells_only[mat]
                             for mat in mesh.GetMaterials()])
    f += GetIsotropicCompressionLinearContribution(fes, E_cells_only, nu, isotropic_cell_contraction)
    vertical_symmetry_plane_penalty_y_dict = {"origin_PDMS_yz": 0,
                      "origin_PDMS_xz": 0,
                      "interior_1_between_PDMS_and_cell_bottom_yz": 0,
                      "ridge_PDMS_yz": 0,
                      "cells_bottom_free_vertical_xz": 0,
                      "inner_2_between_PDMS_and_cell_bottomx_xy": 0,
                      "PDMS_bottom_free_surface_xz": 0,
                      "cells_bottom_free_xy": 0,
                      "cells_bottom_symmetry_yz": 0,
                      "cells_bottom_mid_channel_xz": 0,
                      "cells_top_far_end_yz": 0,
                      "cells_top_y_end_xz": 1e20,
                      "cells_top_free_xy": 0,
                      "cells_top_far_end_to_ridge_xz": 0,
                      "cells_top_ridge_end_yz": 0,
                      "cells_top_ridge_end_xz": 0,
                      "cells_top_symmetry_yz": 0,
                      "PDMS_far_end_yz": 0,
                      "PDMS_y_end_xz": 1e20,
                      "inner_3_between_PDMS_and_cell_top_xy": 0,
                      "PDMS_bottom_xy": 0,
                      "ridge_PDMS_xz": 0
                      }
    vertical_symmetry_plane_penalty_y = CoefficientFunction([vertical_symmetry_plane_penalty_y_dict[bound] for bound in mesh.GetBoundaries()])
    # Impose the y on the long face
    a += vertical_symmetry_plane_penalty_y * u[1] * v[1] * ds # to generally impose the penalty on the long vertical face
    f += vertical_symmetry_plane_penalty_y * y_on_long_edge* v[1] * ds # to shift the intended position to y_on_long_edge
    # For consistency, also add imposed global compression on these constraints
    if orientation_stretch==math.pi/4:
        f += vertical_symmetry_plane_penalty_y * (-chip_compression)*y* v[1] * ds
    else:
        if not orientation_stretch==0:
            f += vertical_symmetry_plane_penalty_y * (-chip_compression * math.sin(orientation_stretch)) * y * v[1] * ds
    f.Assemble()
    u = GridFunction(fes)
    a.Assemble()
    # On thsi Ubuntu installation, using
    # a preconditioner makes things difficult, go directly
    BVP(bf=a, lf=f, gf=u, inverse='sparsecholesky')
    fesp = H1(mesh, order=3, dim=9)  # This is to get the finite element space for the flux (aka, strain or stress)
    strain = GridFunction(fesp)  # Grid function for the stress or strain
    strain.Set(epsilon(u))
    stretch=GridFunction(fesp)
    stretch.Set(epsilon(u)+Id(u.dim))
    fesp = H1(mesh, order=3, dim=9)  # This is to get the finite element space for the flux (aka, strain or stress)
    stress = GridFunction(fesp)  # Grid function for the stress or strain
    stress.Set(sigma(epsilon(u),E,nu))
    return_dict = {
        "fes": fes,
        "u": u,
        "mesh": mesh,
        "strain": strain,
        "stretch": stretch,
        "stress": stress
    }
    return return_dict




def get_simulation_optimize_long_edge(rho=1000, E_PDMS=2e6, nu_PDMS=0.49, E_cells=3e4, nu_cells=0.49, g_acc=9.81,isotropic_cell_contraction=0.05):
    delta_y=0.01
    y_0=get_simulation_with_y_on_vertical_long_face(rho=rho, E_PDMS=E_PDMS, nu_PDMS=nu_PDMS, E_cells=E_cells,
                                                    nu_cells=nu_cells, g_acc=g_acc,chip_compression=0,isotropic_cell_contraction=isotropic_cell_contraction,
                                                    ridge_lateral_partial_adhesion=True,
                                                    orientation_stretch=0,y_on_long_edge=0)
    y_dy=get_simulation_with_y_on_vertical_long_face(rho=rho, E_PDMS=E_PDMS, nu_PDMS=nu_PDMS, E_cells=E_cells,
                                                    nu_cells=nu_cells, g_acc=g_acc,chip_compression=0,isotropic_cell_contraction=isotropic_cell_contraction,
                                                    ridge_lateral_partial_adhesion=True,
                                                    orientation_stretch=0,y_on_long_edge=delta_y)
    # We want the yy-component of the stress to be zero such as to have an equilibrium elongation in the y-direction
    # get this component for the two scenarii, offsetting it for the stress causesd by the cells themselves,
    # which we want to compensate
    e = Id(y_0["u"].dim) * isotropic_cell_contraction
    s = sigma(e, E_cells, nu_cells)
    mesh=y_0["mesh"]
    cell_stress=Integrate(s, mesh, definedon=mesh.Boundaries("cells_top_y_end_xz"))[4]
    stress_y_0=Integrate(y_0["stress"],y_0["mesh"],definedon=y_0["mesh"].Boundaries("cells_top_y_end_xz"))[4]+\
               Integrate(y_0["stress"],y_0["mesh"],definedon=y_0["mesh"].Boundaries("PDMS_y_end_xz"))[4]
    stress_y_dy = Integrate(y_dy["stress"], y_dy["mesh"], definedon=y_dy["mesh"].Boundaries("cells_top_y_end_xz"))[4] + \
                 Integrate(y_dy["stress"], y_dy["mesh"], definedon=y_dy["mesh"].Boundaries("PDMS_y_end_xz"))[4]
    # y0 was calculated with y=0, take it from there
    delta_y_final = (stress_y_0+cell_stress)/(stress_y_dy-stress_y_0)*delta_y+0
    print("Long edge y-displacement for 0 yy-stress = ",delta_y_final,"\n")
    return get_simulation_with_y_on_vertical_long_face(rho=rho, E_PDMS=E_PDMS, nu_PDMS=nu_PDMS, E_cells=E_cells,
                                                    nu_cells=nu_cells, g_acc=g_acc,chip_compression=0,isotropic_cell_contraction=isotropic_cell_contraction,
                                                    ridge_lateral_partial_adhesion=True,
                                                    orientation_stretch=0,y_on_long_edge=delta_y_final)








