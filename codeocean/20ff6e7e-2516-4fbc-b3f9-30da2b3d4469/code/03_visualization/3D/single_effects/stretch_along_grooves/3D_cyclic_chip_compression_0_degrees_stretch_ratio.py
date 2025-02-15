
import pyvista as pv

grid = pv.read('/data/evaluation/vtk/3D_cyclic_chip_compression_0.vtk')


grid.plot(off_screen=True,scalars='strain_ratio_xy',
          show_scalar_bar=True, show_axes=False,cpos=[[-5,-1.25,1],[4,1,1],[0,0,1]],cmap='jet',clim=[1,1.5],
          screenshot="/results/supplementary_2/3D/single_effects/stretch_along_grooves/3D_cyclic_chip_compression_0_degrees_stretch_ratio.png")



