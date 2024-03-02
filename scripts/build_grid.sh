#!/bin/bash

metis_exe=/home/artem/projects/metis-5.1.0/build/Linux-x86_64/programs/mpmetis

if [ -f "$metis_exe" ]; then

    if [ -d utils/build/bin ]; then
        echo "utils found!"     

        val=""
        if [ -f "$1".txt ]; then
            echo "generate from empty format"
            val="txt"     
            ./utils/build/bin/netgen_to_metis $1.$val $1 3
        else
            if [ -f "$1".vol ]; then
                echo "generate from mesh format"     
                val="vol"
                ./utils/build/bin/netgen_mesh_to_metis $1.$val $1
            else
                echo "No input file" 
            fi
        fi

        if [ -f $1 ]; then            
            $metis_exe $1 $2     
            rm $1.npart.*
            mv $1.epart.$2 build/metis     
            rm $1
            ./utils/build/bin/renum_grid_by_metis $1.$val build/metis $1.vtk
            ./utils/build/bin/make_all_geo_data build/ $1.vtk
            ./utils/build/bin/get_average_volume $1.vtk
        else
            echo "Fail netgen to metis convert"
        fi
    else
        echo "utils not found!"     
    fi
else
    echo "metis not found!"     
fi
