import os
os.system('Xvfb :99 -screen 0 1024x768x24 &')
os.environ['DISPLAY'] = ':99'

import pyvista as pv

grid = pv.read('/data/evaluation/vtk/3D_cyclic_chip_compression_90.vtk')


grid.plot(off_screen=True,scalars='cell_orientation_angle',
          show_scalar_bar=True, show_axes=False,cpos=[[-5,-1.25,1],[4,1,1],[0,0,1]],cmap='jet',clim=[0,90],
          screenshot="/results/supplementary_2/3D/single_effects/stretch_perpendicular_to_grooves/3D_cyclic_chip_compression_90_degrees_cell_orientation.png")



