#!/bin/bash

#Построенние картинной плоскости в реальном времени, без кэшировния трассировки
if [ -d utils/build/bin ]; then
    echo "utils found!"         
    
    if [ -d build/bin ]; then        
        #./utils/build/bin/make_griddir_to_plan_projection build/ 2.45 0 0 2 0 0 0.45 0.45 100 100 #для сферы r=1
        ./utils/build/bin/make_griddir_to_plan_projection build/ 3.3 0 0 2 0 0 0.45 0.45 25 25 #для сферы r=0.5 (наблюдатель дальше)
        ./build/bin/get_observer_projection config/directories_cfg.json 
        mpirun -n $1 ./build/bin/tracing_by_observer config/directories_cfg.json
        ./utils/build/bin/rebuild_image_plane build/trace/Solve 1 0.45 0.45 25 25
    else
        echo "build not found"
    fi
else
    echo "utils not found!"     
fi