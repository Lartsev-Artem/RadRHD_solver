#!/bin/bash
cd utils
if [ ! -d build ]; then
    echo "utils not found!"
    cmake -B build -DINSTALL=ON
    cmake --build build -j$1
fi
cd ..
make clean
make prebuild -j$1
make all -j$1
make all -j$1 #первый раз выдаст ошибку из-за параллельной сборки
if [ -d ../grid_data ]; then
    echo "grid_data found!"    
    ./utils/build/bin/make_all_geo_data build/ ../grid_data/sphere_hole23k.vtk
    ./utils/build/bin/get_sphere_direction ../grid_data/directions/surface_sphere_28.vtk build/
fi
