#!/bin/bash

if [ -d utils/build/bin ]; then
    echo "utils found!"
     ./utils/build/bin/netgen_mesh_to_vtk $1.vol $1.vtk
     ./utils/build/bin/netgen_mesh_to_metis $1.vol $1
     /home/artem/projects/metis-5.1.0/build/Linux-x86_64/programs/mpmetis $1 2
     /home/artem/projects/metis-5.1.0/build/Linux-x86_64/programs/mpmetis $1 3
     /home/artem/projects/metis-5.1.0/build/Linux-x86_64/programs/mpmetis $1 4
     rm $1.npart.*
     mv $1.epart.2 build/metis_2
     mv $1.epart.3 build/metis_3
     mv $1.epart.4 build/metis_4
     rm $1
     ./utils/build/bin/make_all_geo_data build/ $1.vtk    
fi
