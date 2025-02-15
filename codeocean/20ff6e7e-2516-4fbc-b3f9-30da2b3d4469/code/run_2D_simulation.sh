#!/usr/bin/env bash
set -ex

# 2D simulations

mkdir -p /results/simulation_output/
mkdir -p /results/simulation_output/2D
mkdir -p /results/simulation_output/2D/self_condensation
mkdir -p /results/simulation_output/2D/cyclic_chip_compression_0_degrees
mkdir -p /results/simulation_output/2D/cyclic_chip_compression_90_degrees


# Change to the directory with the Python simulation files
cd /code/01_simulation/2D_cell_seeding

# First, self condensation by isotropic cellular contraction
python 2D_self_condensation.py > /results/simulation_output/2D/self_condensation/output_2D_self_condensation.txt
# Copy the results to the results section so that they are saved
cp u.sol /results/simulation_output/2D/self_condensation
cp stretch.sol /results/simulation_output/2D/self_condensation
cp strain.sol /results/simulation_output/2D/self_condensation
cp stress.sol /results/simulation_output/2D/self_condensation
# Also copy over the mesh so that the solutions can be drawn and analyzed on the mesh
cp cells_2D_seeded.vol /results/simulation_output/2D/self_condensation


# Stretch along the grooves, same approach
cd /code/01_simulation/2D_cell_seeding

python 2D_stretch_grooves_0.py > /results/simulation_output/2D/cyclic_chip_compression_0_degrees/output_2D_stretch_grooves_0.txt



cp u.sol /results/simulation_output/2D/cyclic_chip_compression_0_degrees
cp stretch.sol /results/simulation_output/2D/cyclic_chip_compression_0_degrees
cp strain.sol /results/simulation_output/2D/cyclic_chip_compression_0_degrees
cp stress.sol /results/simulation_output/2D/cyclic_chip_compression_0_degrees
cp cells_2D_seeded.vol /results/simulation_output/2D/cyclic_chip_compression_0_degrees

# Stretch perpendicular to the grooves, same approach

cd /code/01_simulation/2D_cell_seeding

python 2D_stretch_grooves_90.py > /results/simulation_output/2D/cyclic_chip_compression_90_degrees/output_2D_stretch_grooves_90.txt

cp u.sol /results/simulation_output/2D/cyclic_chip_compression_90_degrees
cp stretch.sol /results/simulation_output/2D/cyclic_chip_compression_90_degrees
cp strain.sol /results/simulation_output/2D/cyclic_chip_compression_90_degrees
cp stress.sol /results/simulation_output/2D/cyclic_chip_compression_90_degrees
cp cells_2D_seeded.vol /results/simulation_output/2D/cyclic_chip_compression_90_degrees




