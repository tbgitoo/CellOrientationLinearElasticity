from ngsolve import *
from ngsolve.bvp import BVP
import netgen.meshing as nm
import math
from CellOrientationLinearElasticity import *




def get_simulation(rho=1000, E=2e6, nu=0, g_acc=9.81,chip_compression=0.2,isotropic_cell_contraction=0.05,
                   ridge_lateral_partial_adhesion=True,orientation_stretch=0):
    m = nm.Mesh(3)
    m.Load("/code/01_simulation/hydrogel/stretching_chip_mesh.vol")
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
    mu = E / 2 / (1 + nu)
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
     # We use here the penalty dictionnary schemes to fix the relevant coordinates
    fes = H1(mesh, order=4, dim=3)
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
    penalty_x_dict = {"mid_chip_symmetry_plane_yz": 1e20,
                      "gel_general_bottom_xy": 0,
                      "end_gel_x_yz": 1e20,
                      "mid_channel_symmetry_plane_y_xz": 0,
                      "gel_free_top_z_xy": 0,
                      "mid_ridge_symmetry_plane_origin_xz": 0,
                      "ridge_far_end_x_yz": 1e20,
                      "rigde_long_vertical_face_y_xz" :1e5, # Partial adhesion
                      "ridge_top_z_xy":0}
    penalty_y_dict  = {"mid_chip_symmetry_plane_yz": 0,
                       "gel_general_bottom_xy": 0,
                       "end_gel_x_yz": 0,
                       "mid_channel_symmetry_plane_y_xz": 1e20,
                       "gel_free_top_z_xy": 0,
                      "mid_ridge_symmetry_plane_origin_xz": 1e20,
                       "ridge_far_end_x_yz": 0,
                       "rigde_long_vertical_face_y_xz" :1e5, # Partial adhesion
                       "ridge_top_z_xy":0}
    penalty_z_dict = {"mid_chip_symmetry_plane_yz": 0,
                      "gel_general_bottom_xy": 1e20,
                      "end_gel_x_yz": 0,
                      "mid_channel_symmetry_plane_y_xz": 0,
                      "gel_free_top_z_xy": 0,
                      "mid_ridge_symmetry_plane_origin_xz": 0,
                      "ridge_far_end_x_yz": 0,
                      "rigde_long_vertical_face_y_xz": 1e5, # Partial adhesion
                      "ridge_top_z_xy": 1e20}
    if not ridge_lateral_partial_adhesion: # We allow retraction only in pseudostatic conditions due to water displacement
        penalty_y_dict["rigde_long_vertical_face_y_xz"]=1e20
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
    f += GetIsotropicCompressionLinearContribution(fes, E, nu, isotropic_cell_contraction)
    f.Assemble()
    u = GridFunction(fes)
    # Use u here to set the Dirichlet boundary condition displacement values
    bound_coef = CoefficientFunction((0, 0, (-z - 0.9)))
    u.Set(bound_coef, BND)
    c = MultiGridPreconditioner(a,inverse='sparsecholesky')
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







