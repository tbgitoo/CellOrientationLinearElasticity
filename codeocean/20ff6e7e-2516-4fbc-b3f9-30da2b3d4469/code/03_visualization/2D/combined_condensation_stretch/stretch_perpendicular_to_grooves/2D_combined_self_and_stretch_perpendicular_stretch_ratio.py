
import pyvista as pv

grid = pv.read('/data/evaluation/vtk/2D_combined_condensation_stretch_90_degrees.vtk')


grid.plot(off_screen=True,scalars='stretch_ratio',
          show_scalar_bar=True, show_axes=False,cpos=[[-3.5,-1.5,1.25],[1,-0.2,0],[0,0,1]],cmap='jet',clim=[1,1.5],
          screenshot="/results/supplementary_2/2D/combined_condensation_stretch/stretch_perpendicular_to_grooves/2D_combined_condensation_stretch_90_degrees_stretch_ratio.png")

