# rm -r build
# cmake -B build
# cmake --build ./build

cmake_minimum_required(VERSION 3.10) 

set(CMAKE_CXX_COMPILER mpic++)

project(solver VERSION 1.0) 

if(NOT CMAKE_CXX_EXTENSIONS)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

include("cmake_modules.cmake")   

FIND_PACKAGE(OpenMP)

add_definitions(-DBUILD_GRAPH -DUSE_MPI -DDEBUG)

add_compile_options(    
    -fopenmp    
    -g           
    # -Og

    # -Werror
    # -pedantic-errors

    # -Wall
    # -Wextra    
    # -Wpedantic

    # -Wcast-align
    
    # -Wconversion    
    # -Wenum-compare
    # -Wfloat-equal    
    # -Woverloaded-virtual
    # -Wredundant-decls    
    # -Wsign-promo


    # json. \todo: add ignore local diagnostic
    -Wno-unused-function

    # MPI warning
    #-Wcast-qual
    -Wno-cast-function-type 

    # Eigen warning
    # -Wctor-dtor-privacy
    # -Wnon-virtual-dtor   
    #-Wsign-conversion
)

include_directories(${EIGEN_DIR}) # библиотека Eigen


set(MAIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(MAIN_DIR_INC "${MAIN_DIR}/include")
set(MAIN_DIR_SRC "${MAIN_DIR}/src")

set(GRAPH_DIR "${MAIN_DIR}/build_graph")
set(LIB_DIR "${MAIN_DIR}/lib")
#set(UTILS_DIR "${MAIN_DIR}/utils")

set(SOURCES_DIR ${MAIN_DIR_INC} ${MAIN_DIR_SRC} ${GRAPH_DIR} ${LIB_DIR} ) #${UTILS_DIR}

FOREACH(dir_to_inc ${SOURCES_DIR})
    INCLUDE_SRC_IN_PRJ(${dir_to_inc} src_files)
ENDFOREACH()

file( GLOB SOURCES ${src_files})

add_subdirectory(lib)

add_executable(${PROJECT_NAME} ${SOURCES})

if(OpenMP_CXX_FOUND)
  target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()


target_link_libraries(${PROJECT_NAME} PUBLIC lib)
