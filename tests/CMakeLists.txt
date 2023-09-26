cmake_minimum_required(VERSION 3.10) 

set(CMAKE_CXX_COMPILER g++)

project(test VERSION 1.0) 

if(NOT CMAKE_CXX_EXTENSIONS)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

include("../cmake_modules.cmake")   

FIND_PACKAGE(OpenMP)

add_compile_options( -fopenmp)

add_definitions(-DBUILD_GRAPH  -DDEBUG)

set(MAIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/..")
set(MAIN_DIR_INC "${MAIN_DIR}/include")
#set(MAIN_DIR_SRC "${MAIN_DIR}/src")

include_directories(${MAIN_DIR})

set(GRAPH_DIR "${MAIN_DIR}/build_graph")
set(LIB_DIR "${MAIN_DIR}/lib")
set(UTILS_DIR "${MAIN_DIR}/utils")
set(TEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

set(SOURCES_DIR ${MAIN_DIR_INC} ${GRAPH_DIR} ${LIB_DIR} ${UTILS_DIR} ${TEST_DIR})

FOREACH(dir_to_inc ${SOURCES_DIR})
    INCLUDE_SRC_IN_PRJ(${dir_to_inc} src_files)
ENDFOREACH()

file( GLOB SOURCES ${src_files})

add_subdirectory("${MAIN_DIR}/lib" "${MAIN_DIR}/lib")

################ Исполняемые файлы
# set(src_files "")
# INCLUDE_SRC_IN_PRJ(${TEST_DIR} src_files)
# file( GLOB EXECUTE_FILES ${src_files})

# FOREACH(f ${EXECUTE_FILES})
#     get_filename_component(name ${f} NAME_WLE) 
    
#     add_executable(${name} ${SOURCES} ${f})
#     message(WARNING ${name})

#     target_link_libraries(${name} PUBLIC lib)

# if(OpenMP_CXX_FOUND)
#   target_link_libraries(${name} PUBLIC OpenMP::OpenMP_CXX)
# endif()


# ENDFOREACH()
#################
# FOREACH(f ${SOURCES})
# message(WARNING ${f})
# ENDFOREACH()

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} PUBLIC lib)

if(OpenMP_CXX_FOUND)
  target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

