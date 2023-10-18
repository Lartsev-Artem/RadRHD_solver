# defining project config 
#DEFCONF 		= SOLVERS LINUX RHLLC DEBUG #BUILD_GRAPH MAKE_TRACE ILLUM USE_MPI  LINUX USE_CUDA
DEFCONF 		= SOLVERS DEBUG BUILD_GRAPH MAKE_TRACE ILLUM USE_MPI LINUX USE_CUDA

# defining working directories
LIB_DIR = lib lib/json lib/files_sys/include lib/Eigen lib/geometry/include lib/mpi_extension
LIB_SRC = lib/files_sys/src lib/geometry lib/mpi_extension lib/json  lib/geometry/src

CUDA_INCDIR = cuda/include cuda/include/interface cuda/include/short_characteristics cuda/include/ray_tracing
CUDA_SRCDIR = cuda/src cuda/src/interface cuda/src/short_characteristics cuda/src/ray_tracing

SOLVERS_DIR = solvers solvers/illum/include solvers/rhllc/include solvers/ray_tracing/include solvers/illum/include/mpi
SOLVERS_SRC = solvers solvers/illum/src solvers/rhllc/src solvers/ray_tracing/src solvers/illum/src/mpi

SRCDIR          = src graph/src make_trace/src ${LIB_SRC} ${SOLVERS_SRC}
INCLUDESDIR     = include graph/include make_trace/include  ${LIB_DIR}  ${CUDA_INCDIR} ${SOLVERS_DIR}
BUILDDIR        = build
OBJDIR          = $(BUILDDIR)/objs
DEPDIR          = $(BUILDDIR)/dep

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
#######################################################################
################ CONFIGURING THE COMPILER #############################
#######################################################################

CXX             = mpic++
CPPFLAGS        = $(DEF_SET) -fopenmp  -Ofast -fPIE
CXXFLAGS        = -std=c++17 #-g #-Wall -Wextra -std=c++11

NVCC 				= nvcc
NVCC_OTHER_FLAGS 	= -Xcompiler "-fopenmp" --expt-relaxed-constexpr #-gencode arch=compute_70,code=sm_70 #совместимость с Eigen3
NVCC_FLAGS 			= $(DEF_SET) -O2  -dc $(NVCC_OTHER_FLAGS)  #-gencode arch=compute_52,code=sm_52

PROGRAM         = run
#TARGET_ARCH	=

COMPILE.cpp     = $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDE_DIRS) $(TARGET_ARCH)
COMPILE.cu     	= $(NVCC) $(NVCC_FLAGS) $(NVCC_OTHER_FLAGS) $(INCLUDE_DIRS) $(TARGET_ARCH)
LINK.cpp 		= $(NVCC) -ccbin=$(CXX) $(NVCC_OTHER_FLAGS)

all: $(PROGRAM)


#Make executable file. Early was $(addprefix $(OBJDIR)/, $^)
$(PROGRAM): %: $(OBJS)
	$(LINK.cpp) $(INCLUDE_DIRS) $(addprefix ./, $^) $(LOADLIBES) $(LDLIBS) -o $@ 
	mv $(PROGRAM) $(BUILDDIR)
	mkdir -p $(BUILDDIR)/graph
	mkdir -p $(BUILDDIR)/illum_geo
	mkdir -p $(BUILDDIR)/Solve
	mkdir -p $(BUILDDIR)/trace

# Building rule for .o files and its .c/.cpp in combination with all .h
%.o: %.cpp
	$(COMPILE.cpp) -c $<
	mv $@ $(OBJDIR)/$@

# Creates the dependecy rules
%.d: %.cpp
	mkdir -p $(OBJDIR)
	mkdir -p $(DEPDIR)
	$(COMPILE.cpp) $^ -MM -MT $(addprefix $(OBJDIR)/, $(@:.d=.o)) > $(DEPDIR)/$@

%.o: %.cu
	$(COMPILE.cu) -c $<
	mv $@ $(OBJDIR)/$@

%.d: %.cu
	mkdir -p $(OBJDIR)
	mkdir -p $(DEPDIR)
	$(COMPILE.cu) $^ -MM -MT $(addprefix $(OBJDIR)/, $(@:.d=.o)) > $(DEPDIR)/$@	

# Includes all .h files
include $(DEP_FILES)
include $(CUDA_DEP_FILES)

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
## попытка собирать только утилиты. с линковкой к библиотеке
# SRCDIR_UTILS         = src lib/files_sys/src lib/geometry lib/mpi_extension lib/json lib/geometry/src utils/netgen_to_vtk
# SRCS_UTILS 			= $(foreach dir,$(SRCDIR_UTILS),$(wildcard $(dir)/*.cpp))
# OBJS_UTILS            = $(patsubst %.cpp, %.o, $(notdir $(SRCS_UTILS)))

# UTIL = util

# utils: $(UTIL)

# $(UTIL): %: $(OBJS_UTILS)
# 	$(LINK.cpp) $(INCLUDE_DIRS) $(addprefix ./, $^) $(LOADLIBES) $(LDLIBS) -o $@ 

# test:
# 	make clean