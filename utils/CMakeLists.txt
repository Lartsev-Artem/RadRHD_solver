cmake_minimum_required(VERSION 3.10) 

project(utils VERSION 1.0) 

get_directory_property(IS_SUBPROJECT PARENT_DIRECTORY)

option(INSTALL "install executable files" OFF)

# set(CMAKE_CXX_COMPILER g++) #компилятор

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(MAIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/..")
set(MAIN_SUB_DIRS "include" "lib")

include_directories("${MAIN_DIR}/include") #путь к файлу prj_config.h
include("${MAIN_DIR}/cmake_modules.cmake")  # функции cmake    


# задаём собственную конфигурацию
if(INSTALL)
    # set(CMAKE_CXX_COMPILER g++) #компилятор
    message(WARNING "installing utils")                   
    add_compile_options(-Ofast -fopenmp -std=c++17)
    add_definitions(-DGENERATE_UTILS)
endif() #INSTALL

#отдельно от общих файлов для генерации набора исполняемых файлов
get_all_files("${CMAKE_CURRENT_SOURCE_DIR}/*.cpp" ${PROJECT_NAME}_SRC_CXX)
get_all_files("${CMAKE_CURRENT_SOURCE_DIR}/*.h" SRC_H)


set (SRC_CXX "")
foreach(sub_dir ${MAIN_SUB_DIRS})	
	get_all_files("${MAIN_DIR}/${sub_dir}/*.cpp" ${PROJECT_NAME}_FILES_CXX)
	get_all_files("${MAIN_DIR}/${sub_dir}/*.h" ${PROJECT_NAME}_FILES_H)
	list(APPEND SRC_CXX ${${PROJECT_NAME}_FILES_CXX})
	list(APPEND SRC_H ${${PROJECT_NAME}_FILES_H})  
endforeach()

# foreach(f ${SRC_H})	
#     message(WARNING ${f})
# endforeach()
get_dirs_from_files("${SRC_H}" ${PROJECT_NAME}_DIR_H)
include_directories(${${PROJECT_NAME}_DIR_H})

#если это подпроект добавить пути в глобальную область видимости
if(IS_SUBPROJECT)
    set(${PROJECT_NAME}_INCLUDE_DIR ${${PROJECT_NAME}_DIR_H} PARENT_SCOPE)
endif()


#подключение библиотек
include_directories(${EIGEN_DIR}) # библиотека Eigen

# if(NOT IS_SUBPROJECT)
#     add_subdirectory("${MAIN_DIR}/lib" "${MAIN_DIR}/lib")    
# endif()
#     include_directories(${lib_INCLUDE_DIR}) # добавляем пути из подпроекта

find_package(VTK REQUIRED)

if(VTK_FOUND)
    add_definitions(-DUSE_VTK) 
else()
    message(WARNING "vtk wasn't found")
endif()

if(INSTALL) # формеруем набор подпрограмм
    foreach(sub_prog ${${PROJECT_NAME}_SRC_CXX})
        get_filename_component(name ${sub_prog}  NAME_WLE) #имя файла без пути и без разширения
        add_executable(${name} ${sub_prog} ${SRC_CXX} ${SRC_H})
        # target_link_libraries(${name} PRIVATE lib)
        if(VTK_FOUND)
            target_link_libraries(${name} PRIVATE  ${VTK_LIBRARIES})
        endif()
    endforeach()
else() # создаем библиотеку
    add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SRC_CXX} ${SRC_CXX} ${SRC_H})
    if(VTK_FOUND)
         target_link_libraries(${PROJECT_NAME} PRIVATE  ${VTK_LIBRARIES})
    endif()
    target_link_libraries(${PROJECT_NAME} PRIVATE lib)
endif()


