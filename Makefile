# defining project config 
DEFCONF 		= SOLVERS DEBUG  BUILD_GRAPH MAKE_TRACE ILLUM USE_CUDA USE_MPI RHLLC RAD_RHD
#DEFCONF 		= SOLVERS DEBUG RHLLC #BUILD_GRAPH MAKE_TRACE ILLUM USE_MPI USE_CUDA

# defining working directories
LIB_DIR = lib lib/json lib/files_sys/include lib/Eigen lib/geometry/include lib/mpi_extension lib/aligned_vec lib/physics/include lib/output_interface
LIB_SRC = lib/files_sys/src lib/geometry lib/mpi_extension lib/json  lib/geometry/src lib/physics/src

CUDA_INCDIR = cuda/include cuda/include/interface cuda/include/short_characteristics cuda/include/ray_tracing cuda/include/short_characteristics/separate_gpu
export CUDA_SRCDIR = cuda/src cuda/src/interface cuda/src/short_characteristics cuda/src/ray_tracing cuda/src/short_characteristics/separate_gpu cuda/src/short_characteristics/separate_gpu/spectrum

SOLVERS_DIR = solvers solvers/illum/include solvers/rhllc/include solvers/ray_tracing/include solvers/illum/include/mpi solvers/illum/include/add_directions solvers/RadRHD/include solvers/illum/include/spectrum
SOLVERS_SRC = solvers solvers/illum/src solvers/rhllc/src solvers/ray_tracing/src solvers/illum/src/mpi solvers/illum/src/add_directions solvers/RadRHD/src solvers/illum/src/mpi/multi_gpu solvers/illum/src/spectrum solvers/illum/src/spectrum/full_spectrum

RESOURCES_DIR = resources

export SRCDIR          = src graph/src make_trace/src ${LIB_SRC} ${SOLVERS_SRC}
INCLUDESDIR     = include graph/include make_trace/include  ${LIB_DIR}  ${CUDA_INCDIR} ${SOLVERS_DIR} build/resources

export ROOT_DIR 	:= $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
export BUILDDIR     = build
export OBJDIR       = $(BUILDDIR)/objs
export DEPDIR       = $(BUILDDIR)/dep
EXEDIR          	= $(BUILDDIR)/bin

SRCEXE 				= src

# specify the list of directories in which the search should be performed.
vpath %.cpp 		$(SRCDIR) 
vpath %.h %.hpp 	$(INCLUDESDIR)
vpath %.o 			$(OBJDIR)
vpath %.d 			$(DEPDIR)
vpath %.cu 			$(CUDA_SRCDIR) 

.SUFFIXES:						# Delete the default suffixes
.SUFFIXES: .cpp .h .o .d .hpp .cu	# Define our suffix list


SRCS 			= $(foreach dir,$(SRCDIR),$(wildcard $(dir)/*.cpp))
OBJS            = $(patsubst %.cpp, %.o, $(notdir $(SRCS)))
DEP_FILES       = $(patsubst %.o, %.d, $(OBJS))

CUDA_SRCS 		= $(foreach dir,$(CUDA_SRCDIR),$(wildcard $(dir)/*.cu))
CUDA_OBJS       = $(patsubst %.cu, %.o, $(notdir $(CUDA_SRCS)))
CUDA_DEP_FILES  = $(patsubst %.o, %.d, $(CUDA_OBJS))

OBJS            += $(patsubst %.cu, %.o, $(notdir $(CUDA_SRCS))) # add cuda objs

INCLUDE_DIRS    = $(addprefix -I ,$(INCLUDESDIR))
DEF_SET 		= $(addprefix -D , $(DEFCONF))

EXE_SRCS 		= $(foreach dir,$(SRCEXE),$(wildcard $(dir)/*.cpp))
EXE_OBJ	 		= $(patsubst %.cpp, %.o, $(notdir $(EXE_SRCS)))
EXE		 		= $(patsubst %.cpp, %.exe, $(notdir $(EXE_SRCS)))
LINK_SRC		= $(filter-out $(EXE_OBJ),$(OBJS))
#######################################################################
################ CONFIGURING THE COMPILER #############################
#######################################################################

# icpc:  -fast -O3 -xHost -ipo
# gcc:  -Ofast march=cpu-type    -flto (-fwhole-program)
export CXX      = mpic++
CPPFLAGS        = $(DEF_SET) -fopenmp -fPIE -Ofast
CXXFLAGS        = -std=c++17 #-g #-Wall -Wextra -std=c++11

export NVCC 		= nvcc
NVCC_OTHER_FLAGS 	= -Xcompiler "-fopenmp" --expt-relaxed-constexpr #-gencode arch=compute_70,code=sm_70 #совместимость с Eigen3
NVCC_FLAGS 			= $(DEF_SET) -O2  -dc $(NVCC_OTHER_FLAGS)  #-gencode arch=compute_52,code=sm_52

PROGRAM         = run
#TARGET_ARCH	=

export COMPILE.cpp	= $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDE_DIRS) $(TARGET_ARCH)
export COMPILE.cu	= $(NVCC) $(NVCC_FLAGS) $(NVCC_OTHER_FLAGS) $(INCLUDE_DIRS) $(TARGET_ARCH)
export LINK.cpp 	= $(NVCC) -ccbin=$(CXX) $(NVCC_OTHER_FLAGS)

all: $(PROGRAM)	

#Make executable file. Early was $(addprefix $(OBJDIR)/, $^)
$(PROGRAM): %: $(OBJS) $(EXE)
#	$(LINK.cpp) $(INCLUDE_DIRS) $(addprefix ./, $^) $(LOADLIBES) $(LDLIBS) -o $@ 
#	mv $(PROGRAM) $(BUILDDIR)
	mkdir -p $(BUILDDIR)/graph
	mkdir -p $(BUILDDIR)/illum_geo
	mkdir -p $(BUILDDIR)/Solve
	mkdir -p $(BUILDDIR)/trace
	mkdir -p $(BUILDDIR)/add_dir

# actions before the building process,
prebuild:
	mkdir -p $(OBJDIR)
	mkdir -p $(DEPDIR)
	mkdir -p $(EXEDIR)	
	$(MAKE) -f $(RESOURCES_DIR)/makefile gen_header_value SRC_FILE=$(ROOT_DIR)/lib/global_consts.h		
	$(MAKE) -f make_dependens	 

# Linking executable files files from all .o in combination with all .h
%.exe: %.cpp
	$(LINK.cpp) $(INCLUDE_DIRS) $(addprefix ./$(OBJDIR)/, $(LINK_SRC)) $(patsubst %.exe, ./$(OBJDIR)/%.o,$@) $(LOADLIBES) $(LDLIBS) -o $(patsubst %.exe, %,$@)
	mv $(patsubst %.exe, %,$@) $(EXEDIR)

# Building rule for .o files and its .c/.cpp in combination with all .h
%.o: %.cpp
	$(COMPILE.cpp) -c $<
	mv $@ $(OBJDIR)/$@

%.o: %.cu
	$(COMPILE.cu) -c $<
	mv $@ $(OBJDIR)/$@

# Creates the dependecy rules
#%.d: %.cpp
#	mkdir -p $(OBJDIR)
#	mkdir -p $(DEPDIR)
#	$(COMPILE.cpp) $^ -MM -MT $(addprefix $(OBJDIR)/, $(@:.d=.o)) > $(DEPDIR)/$@

#%.d: %.cu
#	mkdir -p $(OBJDIR)
#	mkdir -p $(DEPDIR)
#	$(COMPILE.cu) $^ -MM -MT $(addprefix $(OBJDIR)/, $(@:.d=.o)) > $(DEPDIR)/$@	

# Includes all .h files
#include $(DEP_FILES)
#include $(CUDA_DEP_FILES)


# Cleans complete project
.PHONY: clean
clean:
	@echo $(RM)
	$(RM) $(BUILDDIR) -r
	$(RM) *~ 
	$(RM) $(PROGRAM)
	$(RM) "File_Logs.txt"

.PHONY: clean_o
clean_o:
	$(RM) $(OBJDIR) -r
	mkdir -p $(OBJDIR)
