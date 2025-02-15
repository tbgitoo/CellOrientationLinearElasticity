#!/usr/bin/env bash
set -ex

mkdir -p /results/supplementary_2

mkdir -p /results/supplementary_2/3D

mkdir -p /results/supplementary_2/3D/single_effects

mkdir -p /results/supplementary_2/3D/single_effects/stretch_along_grooves

mkdir -p /results/supplementary_2/3D/single_effects/stretch_perpendicular_to_grooves/

mkdir -p /results/supplementary_2/3D/single_effects/self_condensation

mkdir -p /results/supplementary_2/3D/combined_condensation_stretch/

mkdir -p /results/supplementary_2/3D/combined_condensation_stretch/stretch_along_grooves

mkdir -p /results/supplementary_2/3D/combined_condensation_stretch/stretch_perpendicular_to_grooves

# The video buffer should be run in a non-blocking mode (background process), and
# so we start it as a child process (with & at the end of the command)
Xvfb :1 -screen 1 320x240x8&
# Get the process id, we will need to kill this independent process once we're done
pid1=$!

sleep 1

# The pyvista plotting waits for a user to close the
# window, which will never happen in a headless environment, so run it non-blocking
DISPLAY=:1 python -u /code/03_visualization/3D/single_effects/stretch_along_grooves/3D_cyclic_chip_compression_0_degrees.py&

pid2=$!


sleep 20

DISPLAY=:1 python -u /code/03_visualization/3D/single_effects/stretch_along_grooves/3D_cyclic_chip_compression_0_degrees_stretch_ratio.py&

pid2=$!


sleep 20




DISPLAY=:1 python -u /code/03_visualization/3D/single_effects/stretch_perpendicular_to_grooves/3D_cyclic_chip_compression_90_degrees.py&

pid2=$!


sleep 20


DISPLAY=:1 python -u /code/03_visualization/3D/single_effects/stretch_perpendicular_to_grooves/3D_cyclic_chip_compression_90_degrees_stretch_ratio.py&




# The plotting needs some time, if we terminate too early,
# we don't get the output
sleep 20


DISPLAY=:1 python -u /code/03_visualization/3D/single_effects/self_condensation/3D_self_condensation_orientation_angle.py&




# The plotting needs some time, if we terminate too early,
# we don't get the output
sleep 20



DISPLAY=:1 python -u /code/03_visualization/3D/single_effects/self_condensation/3D_self_condensation_stretch_ratio.py&




# The plotting needs some time, if we terminate too early,
# we don't get the output
sleep 20




DISPLAY=:1 python -u /code/03_visualization/3D/combined_condensation_stretch/stretch_along_grooves/3D_combined_self_and_stretch_along_cell_orientation.py&




# The plotting needs some time, if we terminate too early,
# we don't get the output
sleep 20

DISPLAY=:1 python -u /code/03_visualization/3D/combined_condensation_stretch/stretch_along_grooves/3D_combined_self_and_stretch_along_stretch_ratio.py&




# The plotting needs some time, if we terminate too early,
# we don't get the output
sleep 20


DISPLAY=:1 python -u /code/03_visualization/3D/combined_condensation_stretch/stretch_perpendicular_to_grooves/3D_combined_self_and_stretch_perpendicular_cell_orientation.py&




# The plotting needs some time, if we terminate too early,
# we don't get the output
sleep 20


DISPLAY=:1 python -u /code/03_visualization/3D/combined_condensation_stretch/stretch_perpendicular_to_grooves/3D_combined_self_and_stretch_perpendicular_stretch_ratio.py&


# The plotting needs some time, if we terminate too early,
# we don't get the output
sleep 20






kill $pid1


