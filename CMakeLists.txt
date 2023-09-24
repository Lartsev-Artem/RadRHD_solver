# rm -r build
# cmake -B build
# cmake --build ./build

cmake_minimum_required(VERSION 3.2) 

set(CMAKE_CXX_COMPILER mpic++)

project(sovler VERSION 1.0) 

add_compile_options(
    -Werror
    -fopenmp
    -g
    # -Wall
    # -Wextra

    
    # -Wpedantic

    # -Wcast-align
    # -Wcast-qual
    # -Wconversion
    # -Wctor-dtor-privacy
    # -Wenum-compare
    # -Wfloat-equal
    # -Wnon-virtual-dtor
    # -Wold-style-cast
    # -Woverloaded-virtual
    # -Wredundant-decls
    # -Wsign-conversion
    # -Wsign-promo
)

# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_GENERATOR_PLATFORM x64)
   
set(VTK_DIR "/home/artem/projects/VTK9.3/VTK-build") 
find_package(VTK REQUIRED) #include(${VTK_USE_FILE}) 

set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(GRAPH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/build_graph/src")
set(LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib")

set(INC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(INC_GRAPH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/build_graph/include")
set(INC_LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib")


include_directories(${INC_DIR})
include_directories(${INC_GRAPH_DIR})
include_directories(${INC_LIB_DIR})

include_directories("${INC_LIB_DIR}/files_sys/include")
include_directories("${INC_LIB_DIR}/geometry")
include_directories("${INC_LIB_DIR}/json")
include_directories("${INC_LIB_DIR}/mpi_extension")
include_directories("${INC_LIB_DIR}/solvers_types")
include_directories("${INC_LIB_DIR}/Eigen")


file( GLOB SOURCES
    "${SRC_DIR}/*.cpp"
    "${GRAPH_DIR}/*.cpp"
    "${LIB_DIR}/files_sys/src/*.cpp"
    "${LIB_DIR}/geometry/*.cpp"
    "${LIB_DIR}/json/*.cpp"
    "${LIB_DIR}/mpi_extension/*.cpp"
)
add_executable(${PROJECT_NAME} ${SOURCES})


target_link_libraries(${PROJECT_NAME}  ${VTK_LIBRARIES})
