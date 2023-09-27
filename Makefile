# defining project config 
DEFCONF 		= SOLVERS BUILD_GRAPH USE_MPI DEBUG LINUX

# defining working directories
SRCDIR          = src lib/files_sys/src lib/geometry lib/mpi_extension lib/json graph/src lib/geometry/src utils
INCLUDESDIR     = include lib utils lib/json lib/files_sys/include lib/Eigen lib/geometry/include lib/mpi_extension lib/solvers_types graph/include
BUILDDIR        = build
OBJDIR          = $(BUILDDIR)/objs

# specify the list of directories in which the search should be performed.
vpath %.cpp 		$(SRCDIR) 
vpath %.h %.hpp 	$(INCLUDESDIR)
vpath %.o 			$(OBJDIR)
vpath %.d 			$(BUILDDIR)

.SUFFIXES:						# Delete the default suffixes
.SUFFIXES: .cpp .h .o .d .hpp	# Define our suffix list


SRCS 			= $(foreach dir,$(SRCDIR),$(wildcard $(dir)/*.cpp))
OBJS            = $(patsubst %.cpp, %.o, $(notdir $(SRCS)))
DEP_FILES       = $(patsubst %.o, %.d, $(OBJS))
INCLUDE_DIRS    = $(addprefix -I ,$(INCLUDESDIR))
DEF_SET 		= $(addprefix -D , $(DEFCONF))

#######################################################################
################ CONFIGURING THE COMPILER #############################
#######################################################################

CXX             = mpic++
CPPFLAGS        = $(DEF_SET) -fopenmp #-fPIE -Ofast -fopenmp
CXXFLAGS        = -g #-Wall -Wextra -std=c++11

NVCC 				= nvcc
NVCC_OTHER_FLAGS 	= -Xcompiler "-fopenmp"
NVCC_FLAGS 			= $(DEF_SET) -O2 -gencode arch=compute_70,code=sm_70 -dc $(NVCC_OTHER_FLAGS)

PROGRAM         = run
#TARGET_ARCH	=

COMPILE.cpp     = $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDE_DIRS) $(TARGET_ARCH)

all: $(PROGRAM)

#Make executable file. Early was $(addprefix $(OBJDIR)/, $^)
$(PROGRAM): %: $(OBJS)
	$(LINK.cpp) $(INCLUDE_DIRS) $(addprefix ./, $^) $(LOADLIBES) $(LDLIBS) -o $@ 


# Building rule for .o files and its .c/.cpp in combination with all .h
%.o: %.cpp
	$(COMPILE.cpp) -c -o $(OBJDIR)/$@ $<

# Creates the dependecy rules
%.d: %.cpp
	mkdir -p $(OBJDIR)
	$(COMPILE.cpp) $^ -MM -MT $(addprefix $(OBJDIR)/, $(@:.d=.o)) > $(BUILDDIR)/$@

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