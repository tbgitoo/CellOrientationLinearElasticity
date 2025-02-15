mkdir -p /results/evaluation

mkdir -p /results/evaluation/3D

mkdir -p /results/evaluation/vtk

mkdir -p /results/evaluation/3D/self_condensation_and_cyclic_0

mkdir -p /results/evaluation/3D/self_condensation_and_cyclic_90

python /code/02_evaluation/3D_evaluate_orientation_single_effects.py > /results/evaluation/3D/output_3D_evaluate_orientation_single_effects.txt



python /code/02_evaluation/3D_evaluate_orientation_combined_0.py > /results/evaluation/3D/output_3D_evaluate_orientation_combined_0.txt

python /code/02_evaluation/3D_evaluate_orientation_combined_90.py > /results/evaluation/3D/output_3D_evaluate_orientation_combined_90.txt

    