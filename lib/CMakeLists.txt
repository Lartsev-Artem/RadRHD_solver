# rm -r build
# cmake -B build
# cmake --build ./build

cmake_minimum_required(VERSION 3.10) 

set(CMAKE_CXX_COMPILER mpic++)

project(lib VERSION 1.0) 

if(NOT CMAKE_CXX_EXTENSIONS)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()


add_compile_options(    
    -fopenmp    
    -g       

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

# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_GENERATOR_PLATFORM x64)

include_directories("${MAIN_DIR}/lib/Eigen")

#set(VTK_DIR "/home/artem/projects/VTK9.3/VTK-build") 
#find_package(VTK REQUIRED) 

set(MAIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(LIB_DIR "${MAIN_DIR}")
set(INC_DIR "${MAIN_DIR}/../include")

include_directories(${INC_DIR})

set(SOURCES_DIR ${MAIN_DIR}  ${INC_DIR} ${LIB_DIR} )


#поиск поддиректорий
MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child} AND NOT "${child}" STREQUAL "Eigen")
      LIST(APPEND dirlist ${child})
      ENDIF()   
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

#включение каталог до 3 уровней
MACRO(INC_SUBDIR  curdir) 
    list(APPEND subdirs_list "${curdir}")        
    SUBDIRLIST(SUBDIRS ${curdir})
    FOREACH(child ${SUBDIRS})        
        list(APPEND subdirs_list "${curdir}/${child}")            
        SUBDIRLIST(SUBDIRS2 ${curdir}/${child})
        FOREACH(child2 ${SUBDIRS2})            
             list(APPEND subdirs_list "${curdir}/${child}/${child2}")    
            SUBDIRLIST(SUBDIRS3 ${curdir}/${child}/${child2})
            FOREACH(child3 ${SUBDIRS3})               
                list(APPEND subdirs_list "${curdir}/${child}/${child2}/${child3}")                    
            ENDFOREACH()
        ENDFOREACH()        
    ENDFOREACH()
ENDMACRO()

# перебор всех каталогов с исходниками
FOREACH(dir_src ${SOURCES_DIR})
    INC_SUBDIR(${dir_src})
ENDFOREACH()

set(src_files "")
FOREACH(dir ${subdirs_list})
    include_directories("${dir}")
    list(APPEND src_files "${dir}/*.cpp")        
    # message(WARNING ${dir})
ENDFOREACH()  

file( GLOB SOURCES ${src_files})

add_library(${PROJECT_NAME} STATIC ${SOURCES})

#target_link_libraries(${PROJECT_NAME}  ${VTK_LIBRARIES})

FIND_PACKAGE(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

##! Не рабоатет

# # Рекурсивная функция для добавления всех подкаталогов
# function(add_subdirs subdir depth exclude_list)
#   # Получаем список всех элементов в текущем каталоге
#   file(GLOB items RELATIVE ${MAIN_DIR}/${subdir} ${MAIN_DIR}/${subdir}/*)
  
#   # Перебираем все элементы
#   foreach(item ${items})
#     if(IS_DIRECTORY ${MAIN_DIR}/${subdir}/${item})
#       # Если это действительно каталог и глубина не превышает заданное значение,
#       # и каталог не находится в списке исключений, то добавляем его к списку подкаталогов 
#       # и рекурсивно вызываем эту же функцию для каждого вложенного каталога      
#       if((depth LESS_EQUAL ${MAX_DEPTH}) AND (NOT item IN_LIST exclude_list))
#         list(APPEND subdirs_list ${subdir}/${item})
#         MATH(EXPR depth "${depth}+1")

#         include_directories("${MAIN_DIR}/${subdir}/${item}")
#         message(WARNING ${MAIN_DIR}/${subdir}/${item} )

#         add_subdirs(${subdir}/${item} ${depth} exclude_list)
#       endif()    
#     endif()
#   endforeach()
# endfunction()

# # Список каталогов, которые необходимо исключить
# set(EXCLUDE_LIST "Eigen")
# set(MAX_DEPTH 3) # задаем максимальную глубину вложенности

# # Вызываем функцию для начального каталога с начальной глубиной 0 и списком исключений
# add_subdirs("lib" 0 ${EXCLUDE_LIST})
# add_subdirs("build_graph" 0 ${EXCLUDE_LIST})
# add_subdirs("include" 0 ${EXCLUDE_LIST})
