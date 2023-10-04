# defining project config 
DEFCONF 		= SOLVERS BUILD_GRAPH MAKE_TRACE ILLUM USE_MPI DEBUG LINUX USE_CUDA

# defining working directories
LIB_DIR = lib lib/json lib/files_sys/include lib/Eigen lib/geometry/include lib/mpi_extension
LIB_SRC = lib/files_sys/src lib/geometry lib/mpi_extension lib/json  lib/geometry/src

CUDA_DIR = cuda/include cuda/include/interface cuda/include/short_characteristics 
CUDA_SRC = cuda/src cuda/src/interface cuda/src/short_characteristics 

SRCDIR          = src graph/src make_trace/src  ${LIB_SRC} solvers/illum/src solvers #${CUDA_SRC}
INCLUDESDIR     = include graph/include make_trace/include  ${LIB_DIR} solvers/illum/include solvers ${CUDA_DIR}
BUILDDIR        = build
OBJDIR          = $(BUILDDIR)/objs
DEPDIR          = $(BUILDDIR)/dep

# specify the list of directories in which the search should be performed.
vpath %.cpp 		$(SRCDIR) 
vpath %.h %.hpp 	$(INCLUDESDIR)
vpath %.o 			$(OBJDIR)
vpath %.d 			$(DEPDIR)
vpath %.cu 			$(CUDA_SRC) 

.SUFFIXES:						# Delete the default suffixes
.SUFFIXES: .cpp .h .o .d .hpp .cu	# Define our suffix list


SRCS 			= $(foreach dir,$(SRCDIR),$(wildcard $(dir)/*.cpp))
OBJS            = $(patsubst %.cpp, %.o, $(notdir $(SRCS)))
DEP_FILES       = $(patsubst %.o, %.d, $(OBJS))
INCLUDE_DIRS    = $(addprefix -I ,$(INCLUDESDIR))
DEF_SET 		= $(addprefix -D , $(DEFCONF))

#######################################################################
################ CONFIGURING THE COMPILER #############################
#######################################################################

CXX             = mpic++
CPPFLAGS        = $(DEF_SET) -fopenmp  -Ofast #-fPIE
CXXFLAGS        = #-g #-Wall -Wextra -std=c++11

NVCC 				= nvcc
NVCC_OTHER_FLAGS 	= -Xcompiler "-fopenmp" --expt-relaxed-constexpr #совместимость с Eigen3
NVCC_FLAGS 			= $(DEF_SET) -O2 -gencode arch=compute_70,code=sm_70 -dc $(NVCC_OTHER_FLAGS)

PROGRAM         = run
#TARGET_ARCH	=

COMPILE.cpp     = $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDE_DIRS) $(TARGET_ARCH)


all: $(PROGRAM)


#Make executable file. Early was $(addprefix $(OBJDIR)/, $^)
$(PROGRAM): %: $(OBJS)
	$(LINK.cpp) $(INCLUDE_DIRS) $(addprefix ./, $^) $(LOADLIBES) $(LDLIBS) -o $@ 
	mv $(PROGRAM) $(BUILDDIR)

# Building rule for .o files and its .c/.cpp in combination with all .h
%.o: %.cpp
	$(COMPILE.cpp) -c -o $(OBJDIR)/$@ $<

# Creates the dependecy rules
%.d: %.cpp #%.cu
	mkdir -p $(OBJDIR)
	mkdir -p $(DEPDIR)
	$(COMPILE.cpp) $^ -MM -MT $(addprefix $(OBJDIR)/, $(@:.d=.o)) > $(DEPDIR)/$@

# Includes all .h files
include $(DEP_FILES)

# Cleans complete project
.PHONY: clean
clean:
	@echo $(RM)
	$(RM) $(BUILDDIR) -r
	$(RM) *~ 
	$(RM) $(PROGRAM)
	$(RM) "File_Logs.txt"


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

CUDA_SRCS 			= $(foreach dir,$(CUDA_SRC),$(wildcard $(dir)/*.cu))
CUDA_OBJS            = $(patsubst %.cu, %.o, $(notdir $(CUDA_SRCS)))
CUDA_DEP_FILES       = $(patsubst %.o, %.d, $(CUDA_OBJS))

COMPILE.cu     = $(NVCC) $(NVCC_FLAGS) $(NVCC_OTHER_FLAGS) $(INCLUDE_DIRS) $(TARGET_ARCH)

%.o: %.cu
	$(COMPILE.cu) -c -o $(OBJDIR)/$@ $<

%.d: %.cu
	mkdir -p $(OBJDIR)
	mkdir -p $(DEPDIR)
	$(COMPILE.cu) $^ -MM -MT $(addprefix $(OBJDIR)/, $(@:.d=.o)) > $(DEPDIR)/$@

cuda: $(CUDA_OBJS)
	$(LINK.cu) $(INCLUDE_DIRS) $(addprefix ./, $^) $(LOADLIBES) $(LDLIBS) -o $@ 