#!/bin/bash

#Построенние картинной плоскости c предварительным кэшированием трассировки
if [ -d utils/build/bin ]; then
    echo "utils found!"         
    
    if [ -d build/bin ]; then        
        ./utils/build/bin/make_griddir_to_plan_projection build/ 2.45 0 0 2 0 0 0.45 0.45 25 25        
        mpirun -n $1 ./build/bin/graph_builder config/directories_cfg.json
        mpirun -n $1 ./build/bin/trace_builder config/directories_cfg.json 
        #./build/bin/illum_face_solver config/directories_cfg.json 
        mpirun -n $1 ./build/bin/illum_multi_gpu_solver config/directories_cfg.json 
        ./build/bin/get_observer_projection config/directories_cfg.json 
        ./utils/build/bin/make_observer_plane build/Solve/Solve0Illum build/trace/ray_trace0.bin build/plane_cfg.bin build/trace/plane.vtk
    else
        echo "build not found"
    fi    
else
    echo "utils not found!"     
fi
