cmake_minimum_required(VERSION 3.10) 

project(lib VERSION 1.0) 

get_directory_property(IS_SUBPROJECT PARENT_DIRECTORY)

option(DEBUG_COMPILER "compilation with debugging" OFF)

# задаём собственную конфигурацию
if(NOT IS_SUBPROJECT)
    message(WARNING "building the library separately from the project")
    
    set(CMAKE_CXX_COMPILER g++) #компилятор
        
    include("${CMAKE_CURRENT_SOURCE_DIR}/../cmake_modules.cmake")  # функции cmake    
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../include") #путь к файлу prj_config.h

    if(DEBUG_COMPILER) 
    add_compile_options(   
    -g           
    
    -Werror
    -pedantic-errors

    -Wall
    -Wextra    
    -Wpedantic

    -Wcast-align
    
    -Wconversion    
    -Wenum-compare
    -Wfloat-equal    
    -Woverloaded-virtual
    -Wredundant-decls    
    -Wsign-promo

    # json. \todo: add ignore local diagnostic
    -Wno-unused-function
    
    -Wno-cast-function-type # MPI warning

    # Eigen warning
        # -Wctor-dtor-privacy
        # -Wnon-virtual-dtor   
        # -Wsign-conversion
    )
    else()
        add_compile_options(-Ofast -fopenmp -std=c++17)
    endif()
    
endif() #IS_SUBPROJECT

include_directories(${EIGEN_DIR}) # библиотека Eigen

find_package(VTK REQUIRED)

if(VTK_FOUND)
    add_definitions(-DUSE_VTK) 
else()
    message(WARNING "vtk wasn't found")
endif()

set(SOURCES_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

FOREACH(dir_to_inc ${SOURCES_DIR})
    INCLUDE_SRC_IN_PRJ(${dir_to_inc} src_files)
ENDFOREACH()

file( GLOB SOURCES ${src_files})

add_library(${PROJECT_NAME} STATIC ${SOURCES})

if(VTK_FOUND)
     target_link_libraries(${PROJECT_NAME} PUBLIC  ${VTK_LIBRARIES})
endif()

