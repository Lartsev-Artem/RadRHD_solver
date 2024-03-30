#!/bin/bash

if [ -d utils/build/bin ]; then
    echo "utils found!"         
    
    if [ -d build/bin ]; then
        if [ -f $1 -a -f $2"0.bin" ]; then                    
            
            for (( i=0; i < $3; i++ ))
                do
                     ./utils/build/bin/add_bin_scalar_to_vtk $1 $2$i.bin double illum_$i                    
                     echo add $i
                done
        else
            echo "file $1 OR $2 not found"
        fi
    else
        echo "build not found"
    fi    
else
    echo "utils not found!"     
fi
