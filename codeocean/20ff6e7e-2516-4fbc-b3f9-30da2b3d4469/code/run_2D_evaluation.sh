mkdir -p /results/evaluation
mkdir -p /results/evaluation/2D

mkdir -p /results/evaluation/vtk

mkdir -p /results/evaluation/2D/self_condensation_and_cyclic_0

mkdir -p /results/evaluation/2D/self_condensation_and_cyclic_90



python /code/02_evaluation/2D_evaluate_orientation_single_effects.py > /results/evaluation/2D/output_2D_evaluate_orientation_single_effects.txt

python /code/02_evaluation/2D_evaluate_orientation_combined.py > /results/evaluation/2D/output_2D_evaluate_orientation_combined.txt


    