
import pyvista as pv

grid = pv.read('/data/evaluation/vtk/3D_self_condensation.vtk')


grid.plot(off_screen=True,scalars='strain_ratio_xy',
          show_scalar_bar=True, show_axes=False,cpos=[[-3.5,-1.5,1.25],[1,-0.2,0],[0,0,1]],cmap='jet',clim=[1,1.5],
          screenshot="/results/supplementary_2/3D/single_effects/self_condensation/3D_self_condensation_stretch_ratio.png")