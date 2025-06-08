#!/bin/bash

if [ -d utils/build/bin ]; then
    echo "utils found!"         
    
    if [ -d build/tests ]; then
        if [ -f $1  ]; then
            names=(direct_refl inverse_refl avg_points_inverse_refl lin_inverse_refl avg_cells_refl)
            ./utils/build/bin/get_sphere_direction $1 build/ 1
            ./build/tests/illum/gen_illum_by_dirs config/directories_cfg.json build/illum.bin 0
            ./utils/build/bin/set_bin_scalar_to_vtk $1 build/illum.bin build/new_grid_dir.vtk double src_illum
            for (( i=0; i < 5; i++ ))
                do
                    ./build/tests/illum/bound_reflect_test config/directories_cfg.json build/illum.bin build/ref$i.bin $i                    
                    ./utils/build/bin/add_bin_scalar_to_vtk build/new_grid_dir.vtk build/ref$i.bin double ${names[$i]}
                    rm build/ref$i.bin
                done            
            rm build/illum.bin
        else
            echo "file $1 not found"
        fi
    else
        echo "build not found"
    fi    
else
    echo "utils not found!"     
fi
