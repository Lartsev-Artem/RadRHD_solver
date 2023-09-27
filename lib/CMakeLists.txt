cmake_minimum_required(VERSION 3.10) 

project(lib VERSION 1.0) 

get_directory_property(IS_SUBPROJECT PARENT_DIRECTORY)

option(DEBUG_COMPILER "compilation with debugging" OFF)

include("../cmake_modules.cmake")  # функции cmake    
include_directories("../include") #путь к файлу prj_config.h

# задаём собственную конфигурацию
if(NOT IS_SUBPROJECT)
    message(WARNING "building the library separately from the project")
    
    set(CMAKE_CXX_COMPILER /usr/bin/g++) #компилятор
            
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


get_all_files("${CMAKE_CURRENT_SOURCE_DIR}/*.cpp" ${PROJECT_NAME}_FILES_CXX)
get_all_files("${CMAKE_CURRENT_SOURCE_DIR}/*.h" ${PROJECT_NAME}_FILES_H)

get_dirs_from_files("${${PROJECT_NAME}_FILES_H}" ${PROJECT_NAME}_DIR_H)
include_directories(${PROJECT_NAME} ${${PROJECT_NAME}_DIR_H})

#если это подпроект добавить пути в глобальную область видимости
if(IS_SUBPROJECT)
    set(${PROJECT_NAME}_INCLUDE_DIR ${${PROJECT_NAME}_DIR_H} PARENT_SCOPE)
endif()


#подключение библиотек
include_directories(${EIGEN_DIR}) # библиотека Eigen

find_package(VTK REQUIRED)

if(VTK_FOUND)
    add_definitions(-DUSE_VTK) 
else()
    message(WARNING "vtk wasn't found")
endif()

add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_FILES_CXX}  ${${PROJECT_NAME}_FILES_H})


if(VTK_FOUND)
     target_link_libraries(${PROJECT_NAME} PUBLIC  ${VTK_LIBRARIES})
endif()

