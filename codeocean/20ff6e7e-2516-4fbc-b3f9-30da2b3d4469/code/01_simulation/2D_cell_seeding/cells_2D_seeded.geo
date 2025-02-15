algebraic3d

solid PDMS_base_x1 = plane (0, 0, 0; -1, 0, 0);

solid PDMS_base_y1 = plane (0, 0, 0; 0, -1, 0);

solid PDMS_base_z1 = plane (0, 0, -0.35; 0, 0, -1);

solid PDMS_base_x2 = plane (8, 0, 0; 1, 0, 0);

solid PDMS_base_y2 = plane (0, 0.35, 0; 0, 1, 0);

solid PDMS_base_z2 = plane (0, 0, 0.35; 0, 0, 1);


solid px_removal = plane (-0.5, 0, 0; -1, 0, 0);

solid py_removal = plane (0, -0.5, 0; 0, -1, 0);

solid pz_removal = plane (0, 0, 0; 0, 0, -1);

solid px_removal_end = plane (5, 0, 0; 1, 0, 0);

solid py_removal_end = plane (0, 0.175, 0; 0, 1, 0);

solid pz_removal_end = plane (0, 0, 0.5; 0, 0, 1);



solid cells_top_x1 = plane (0, 0, 0; -1, 0, 0);

solid cells_top_y1 = plane (0, 0, 0; 0, -1, 0);

solid cells_top_z1 = plane (0, 0, 0.35; 0, 0, -1);

solid cells_top_x2 = plane (8, 0, 0; 1, 0, 0);

solid cells_top_y2 = plane (0, 0.35, 0; 0, 1, 0);

solid cells_top_z2 = plane (0, 0, 0.36; 0, 0, 1);


solid cells_bottom_x1 = plane (0, 0, 0; -1, 0, 0);

solid cells_bottom_y1 = plane (0, 0, 0; 0, -1, 0);

solid cells_bottom_z1 = plane (0, 0, 0; 0, 0, -1);

solid cells_bottom_x2 = plane (5, 0, 0; 1, 0, 0);

solid cells_bottom_y2 = plane (0, 0.165, 0; 0, 1, 0);

solid cells_bottom_z2 = plane (0, 0, 0.01; 0, 0, 1);





solid PDMS_block = PDMS_base_x1 and PDMS_base_y1 and PDMS_base_z1 and PDMS_base_x2 and PDMS_base_y2 and PDMS_base_z2;

solid cells_top_layer_full = cells_top_x1 and cells_top_y1 and cells_top_z1 and cells_top_x2 and cells_top_y2 and cells_top_z2;

solid removal_block = px_removal and py_removal and pz_removal and px_removal_end and py_removal_end and pz_removal_end;

solid PDMS_ridge = PDMS_block and not removal_block;

solid cells_top = cells_top_layer_full and not removal_block;

solid cells_bottom = cells_bottom_x1 and cells_bottom_y1 and cells_bottom_z1 and cells_bottom_x2 and cells_bottom_y2 and cells_bottom_z2;

tlo PDMS_ridge;

tlo cells_top;

tlo cells_bottom;
