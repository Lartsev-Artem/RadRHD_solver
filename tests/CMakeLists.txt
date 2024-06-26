cmake_minimum_required(VERSION 3.20) 

set(CMAKE_CXX_COMPILER g++)

project(tests VERSION 1.0) 

if(NOT CMAKE_CXX_EXTENSIONS)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include("../cmake_modules.cmake")   

find_package(OpenMP)

add_compile_options(-fopenmp)

add_definitions(-DBUILD_GRAPH  -DDEBUG -DMAKE_TRACE -DSOLVERS -DILLUM)

set(MAIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/..")
set(MAIN_SUB_DIRS "graph" "include" "lib" "utils" "solvers" "make_trace")

#отдельно от общих файлов для генерации набора исполняемых файлов
get_all_files("${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp" ${PROJECT_NAME}_SRC_CXX)

set (SRC_H "") # включаемые файлы тестов в общий список
set (SRC_CXX "")
foreach(sub_dir ${MAIN_SUB_DIRS})	
	get_all_files("${MAIN_DIR}/${sub_dir}/*.cpp" ${PROJECT_NAME}_FILES_CXX)
	get_all_files("${MAIN_DIR}/${sub_dir}/*.h" ${PROJECT_NAME}_FILES_H)
	list(APPEND SRC_CXX ${${PROJECT_NAME}_FILES_CXX})
	list(APPEND SRC_H ${${PROJECT_NAME}_FILES_H})  
endforeach()

get_dirs_from_files("${SRC_H}" ${PROJECT_NAME}_DIR_H)
include_directories(${${PROJECT_NAME}_DIR_H})

#добавление библиотек
include_directories(${EIGEN_DIR}) # библиотека Eigen

find_package(VTK REQUIRED)

if(VTK_FOUND)
    add_definitions(-DUSE_VTK) 
else()
    message(WARNING "vtk wasn't found")
endif()

# add_subdirectory("${MAIN_DIR}/lib" "${MAIN_DIR}/lib")
# include_directories(${lib_INCLUDE_DIR}) # добавляем пути из подпроекта

# add_subdirectory("${MAIN_DIR}/utils" "${MAIN_DIR}/utils")
# include_directories(${utils_INCLUDE_DIR}) # добавляем пути из подпроекта

# формеруем набор подпрограмм
foreach(sub_prog ${${PROJECT_NAME}_SRC_CXX})
  get_filename_component(name ${sub_prog}  NAME_WLE) #имя файла без пути и без разширения
	add_executable(${name} ${sub_prog} ${SRC_H} ${SRC_CXX})  
  #линковка с библиотеками
  # target_link_libraries(${name} PRIVATE lib)

  # target_link_libraries(${name} PRIVATE utils)
  
  if(OpenMP_CXX_FOUND)
    target_link_libraries(${name} PRIVATE OpenMP::OpenMP_CXX)
  endif()

  if(VTK_FOUND)
     target_link_libraries(${name} PRIVATE  ${VTK_LIBRARIES})
  endif()

endforeach()


