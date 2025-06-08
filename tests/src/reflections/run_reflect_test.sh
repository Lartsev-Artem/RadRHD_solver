#!/bin/bash

if [ -d utils/build/bin ]; then
    echo "utils found!"         
    
    if [ -d build/tests ]; then
        if [ -f $1  ]; then                                                        
            ./utils/build/bin/get_sphere_direction $1 build/ 1
            ./build/tests/illum/gen_illum_by_dirs config/directories_cfg.json build/illum.bin 0
            ./build/tests/illum/bound_reflect_test config/directories_cfg.json build/illum.bin build/ref0.bin 0
            ./build/tests/illum/bound_reflect_test config/directories_cfg.json build/illum.bin build/ref1.bin 1
            ./build/tests/illum/bound_reflect_test config/directories_cfg.json build/illum.bin build/ref2.bin 2
            ./build/tests/illum/bound_reflect_test config/directories_cfg.json build/illum.bin build/ref3.bin 3
            ./build/tests/illum/bound_reflect_test config/directories_cfg.json build/illum.bin build/ref4.bin 4
            ./utils/build/bin/set_bin_scalar_to_vtk $1 build/illum.bin build/new_grid_dir.vtk double src_illum
            ./utils/build/bin/add_bin_scalar_to_vtk build/new_grid_dir.vtk build/ref0.bin double direct_refl
            ./utils/build/bin/add_bin_scalar_to_vtk build/new_grid_dir.vtk build/ref1.bin double inverse_refl
            ./utils/build/bin/add_bin_scalar_to_vtk build/new_grid_dir.vtk build/ref2.bin double avg_points_inverse_refl
            ./utils/build/bin/add_bin_scalar_to_vtk build/new_grid_dir.vtk build/ref3.bin double lin_inverse_refl
            ./utils/build/bin/add_bin_scalar_to_vtk build/new_grid_dir.vtk build/ref4.bin double avg_cells_refl
        else
            echo "file $1 not found"
        fi
    else
        echo "build not found"
    fi    
else
    echo "utils not found!"     
fi
