#!/usr/bin/env bash
set -ex

# First, the 3D hydrogel simulation
# These are the output directories, we simulat
# stretch at 0° and 90° relative to the groove long axis and self condensation
mkdir -p /results/simulation_output/
mkdir -p /results/simulation_output/3D
mkdir -p /results/simulation_output/3D/self_condensation
mkdir -p /results/simulation_output/3D/cyclic_chip_compression_0_degrees
mkdir -p /results/simulation_output/3D/cyclic_chip_compression_90_degrees

# Change to the directory with the Python simulation files
cd /code/01_simulation/hydrogel

# First, self condensation by isotropic cellular contraction
python 3D_self_condensation.py > /results/simulation_output/3D/self_condensation/output_3D_self_condensation.txt

# Copy the results to the results section so that they are saved
cp u.sol /results/simulation_output/3D/self_condensation
cp stretch.sol /results/simulation_output/3D/self_condensation
cp strain.sol /results/simulation_output/3D/self_condensation
cp stress.sol /results/simulation_output/3D/self_condensation
# Also copy over the mesh so that the solutions can be drawn and analyzed on the mesh
cp stretching_chip_mesh.vol /results/simulation_output/3D/self_condensation


# Stretch along the grooves, same approach
python /code/01_simulation/hydrogel/3D_stretch_grooves_0_degrees.py >  /results/simulation_output/3D/cyclic_chip_compression_0_degrees/output_3D_stretch_grooves_0_degrees.txt
cp u.sol /results/simulation_output/3D/cyclic_chip_compression_0_degrees
cp stretch.sol /results/simulation_output/3D/cyclic_chip_compression_0_degrees
cp strain.sol /results/simulation_output/3D/cyclic_chip_compression_0_degrees
cp stress.sol /results/simulation_output/3D/cyclic_chip_compression_0_degrees
cp stretching_chip_mesh.vol /results/simulation_output/3D/cyclic_chip_compression_0_degrees

# Stretch perpendicular to the grooves, same approach
python /code/01_simulation/hydrogel/3D_stretch_grooves_90_degrees.py > /results/simulation_output/3D/cyclic_chip_compression_90_degrees/output_3D_stretch_grooves_90_degrees.txt

cp u.sol /results/simulation_output/3D/cyclic_chip_compression_90_degrees
cp stretch.sol /results/simulation_output/3D/cyclic_chip_compression_90_degrees
cp strain.sol /results/simulation_output/3D/cyclic_chip_compression_90_degrees
cp stress.sol /results/simulation_output/3D/cyclic_chip_compression_90_degrees
cp stretching_chip_mesh.vol /results/simulation_output/3D/cyclic_chip_compression_90_degrees





