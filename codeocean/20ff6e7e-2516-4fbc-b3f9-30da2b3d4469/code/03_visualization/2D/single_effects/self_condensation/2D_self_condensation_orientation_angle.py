
import pyvista as pv

grid = pv.read('/data/evaluation/vtk/2D_self_condensation.vtk')


grid.plot(off_screen=True,scalars='cell_orientation_angle',
          show_scalar_bar=True, show_axes=False,cpos=[[-3.5,-1.5,1.25],[1,-0.2,0],[0,0,1]],cmap='jet',clim=[0,90],
          screenshot="/results/supplementary_2/2D/single_effects/self_condensation/2D_self_condensation_cell_orientation.png")