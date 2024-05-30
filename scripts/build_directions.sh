#!/bin/bash

if [ -d utils/build/bin ]; then
    echo "utils found!"         
    
    if [ -d build/bin ]; then
        if [ -f $1 ]; then                    
            ./utils/build/bin/get_sphere_direction $1 build/
            ./build/bin/graph_builder /home/artem/projects/RadRHD_solver/config/directories_cfg.json
            ./build/bin/trace_builder /home/artem/projects/RadRHD_solver/config/directories_cfg.json
            ./build/bin/get_reflection /home/artem/projects/RadRHD_solver/config/directories_cfg.json
        else
            echo "file sphere direction not found"
        fi
    else
        echo "build not found"
    fi    
else
    echo "utils not found!"     
fi
