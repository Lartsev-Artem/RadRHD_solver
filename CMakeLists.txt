# rm -r build
# cmake -B build
# cmake --build ./build

cmake_minimum_required(VERSION 3.20) 

set(CMAKE_CXX_COMPILER mpic++)

project(solver VERSION 1.0) 

if(NOT CMAKE_CXX_EXTENSIONS)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

include("cmake_modules.cmake")   

set(PRJ_NAME_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(MAIN_SUB_DIRS "include" "src" "lib" "graph" "make_trace" "solvers")

add_definitions(-DBUILD_GRAPH -DUSE_MPI -DRRHD_DEBUG -DMAKE_TRACE -DSOLVERS -DILLUM)

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

#поис всех файлов в подпроектах
set (SRC_CXX "")
set (SRC_H "")
foreach(sub_dir ${MAIN_SUB_DIRS})	
	get_all_files("${PRJ_NAME_DIR}/${sub_dir}/*.cpp" ${PROJECT_NAME}_FILES_CXX)
	get_all_files("${PRJ_NAME_DIR}/${sub_dir}/*.h" ${PROJECT_NAME}_FILES_H)
	list(APPEND SRC_CXX ${${PROJECT_NAME}_FILES_CXX})
	list(APPEND SRC_H ${${PROJECT_NAME}_FILES_H})
endforeach()

#добавление путей в проект
get_dirs_from_files("${SRC_H}" ${PROJECT_NAME}_DIR_H)
include_directories(${${PROJECT_NAME}_DIR_H})

#подключение библиотек
find_package(OpenMP)

include_directories(${EIGEN_DIR}) # библиотека Eigen

# add_subdirectory(lib)
# include_directories(${lib_INCLUDE_DIR}) # добавляем пути из подпроекта

add_executable(${PROJECT_NAME}  ${SRC_CXX}  ${SRC_H})

#линковка с библиотеками
if(OpenMP_CXX_FOUND)
  target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

# target_link_libraries(${PROJECT_NAME} PRIVATE lib)
