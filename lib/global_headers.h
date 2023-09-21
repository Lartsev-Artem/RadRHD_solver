#ifndef GLOBAL_HEADERS
#define GLOBAL_HEADERS

#include <inttypes.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <memory>
#include <string>
#include <vector>

#include "prj_config.h"


#ifdef USE_MPI
#include "mpi.h"
#endif

#ifdef USE_VTK
#include <vtk-9.0/vtkSmartPointer.h>
#include <vtk-9.0\vtkCellArray.h>
#include <vtk-9.0\vtkCellData.h>
#include <vtk-9.0\vtkDataObject.h>
#include <vtk-9.0\vtkDataSet.h>
#include <vtk-9.0\vtkDataSetAttributes.h>
#include <vtk-9.0\vtkDoubleArray.h>
#include <vtk-9.0\vtkGenericDataObjectReader.h>
#include <vtk-9.0\vtkGenericDataObjectWriter.h>
#include <vtk-9.0\vtkIdList.h>
#include <vtk-9.0\vtkLine.h>
#include <vtk-9.0\vtkLineSource.h>
#include <vtk-9.0\vtkMath.h>
#include <vtk-9.0\vtkNamedColors.h>
#include <vtk-9.0\vtkPointData.h>
#include <vtk-9.0\vtkPoints.h>
#include <vtk-9.0\vtkQuad.h>
#include <vtk-9.0\vtkTetra.h>
#include <vtk-9.0\vtkTriangle.h>
#include <vtk-9.0\vtkUnsignedCharArray.h>
#include <vtk-9.0\vtkUnstructuredGrid.h>

#endif  // USE_VTK

#if !defined CLASTER && !defined LINUX
#include <direct.h>
#include <sys/stat.h>

#endif

#if defined BUILD_GRAPH
int RunBuildModule(const std::string& name_file_settings);
#endif  // BUILD_GRAPH

#if defined BUILD_DATA_TO_ILLUM
int RunMakeModule(std::string name_file_settings, int a, int b);
#endif  // BUILD_DATA_TO_ILLUM

#if defined SOLVERS
int RunSolveModule(int argc, char* argv[],
                   const std::string& name_file_settings);
#endif  // SOLVERS

#if defined UTILS
int RunUtilsModule(int argc, char* argv[]);
#endif  // UTILS

#endif  // GLOBAL_HEADERS
