algebraic3d

solid full_sphere = sphere (0, 0, 0; 1);

solid px = plane (0, 0, 0; -1, 0, 0);

solid py = plane (0, 0, 0; 0, -1, 0);

solid pz = plane (0, 0, 0; 0, 0, -1);

solid px_end = plane (8, 0, 0; 1, 0, 0);

solid py_end = plane (0, 0.35, 0; 0, 1, 0);

solid pz_end = plane (0, 0, 2; 0, 0, 1);


solid px_removal = plane (-0.5, 0, 0; -1, 0, 0);

solid py_removal = plane (0, -0.5, 0; 0, -1, 0);

solid pz_removal = plane (0, 0, -0.5; 0, 0, -1);

solid px_removal_end = plane (5, 0, 0; 1, 0, 0);

solid py_removal_end = plane (0, 0.175, 0; 0, 1, 0);

solid pz_removal_end = plane (0, 0, 0.35; 0, 0, 1);




solid removal_block = px_removal and py_removal and pz_removal and px_removal_end and py_removal_end and pz_removal_end;

solid gel_block = px and py and pz and px_end and py_end and pz_end and not removal_block;


tlo gel_block;
