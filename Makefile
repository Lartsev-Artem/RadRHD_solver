## Флаги сборки
KEYS	:= SOLVERS BUILD_GRAPH MAKE_TRACE ILLUM RHLLC RadRHD
KEYS 	+= DEBUG USE_CUDA USE_MPI
DEFINES = $(addprefix -D , $(KEYS))
## Компиляторы и настройки
CONFIG ?= release
CXX 		:= mpiicpx
NVCC 		:= nvcc
CXXFLAGS 	:= $(DEFINES) -fopenmp -fPIE -std=c++17 
NVCCFLAGS 	:= $(DEFINES) --expt-relaxed-constexpr -dc #-gencode arch=compute_70,code=sm_70 #-Xcompiler "-fopenmp"

LDFLAGS 	:= -fopenmp -L/usr/local/cuda/lib64 -lcudart

# Настройки флагов для конфигураций
ifeq ($(CONFIG),debug)
	CXXFLAGS 	+= -g -O0 #-Wall -Wextra -std=c++11
	NVCCFLAGS 	+= -G -g
else
# gcc:  -Ofast march=cpu-type    -flto (-fwhole-program)
	CXXFLAGS 	+= -Ofast #-xHost -ipo
	NVCCFLAGS 	+= -O3
	#LDFLAGS 	+= -ipo
endif

## Пути 
export ROOT_DIR 	:= $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
#export BUILD_DIR = build/$(CONFIG)
export BUILD_DIR 	= build
RESOURCES_DIR 		= resources
BUILD_BIN          	= $(BUILD_DIR)/bin
BUILD_OBJ          	= $(BUILD_DIR)/obj
EXE_DIR 			= src
TESTS_DIR 			= tests/src
TESTS_BIN 			= $(BUILD_DIR)/tests

## Пути к исходным файлам
LIB_SRC_DIR 	= lib/files_sys/src/ lib/geometry/src/ lib/json/src/ lib/mpi_extension/src/ lib/physics/src/ lib/grids/src/
SOLVERS_SRC_DIR = solvers/illum/src/ solvers/hllc/src/ solvers/rhllc/src/ solvers/ray_tracing/src/  solvers/RadRHD/src/

CXX_SRC_DIRS 	= graph/src/ make_trace/src/ ${LIB_SRC_DIR} ${SOLVERS_SRC_DIR}
CUDA_SRC_DIRS 	= cuda/src/

## Пути к заголовочным файлам
LIB_INC_DIR 	= lib lib/aligned_vec lib/files_sys/include lib/geometry/include lib/grids/include lib/json/include lib/mpi_extension/include lib/output_interface lib/physics/include
SOLVERS_INC_DIR = solvers/illum/include/ solvers/hllc/include/ solvers/rhllc/include solvers/ray_tracing/include  solvers/RadRHD/include
CUDALIB_INC_DIR = cuda/include

CUDA_INC_DIRS 	= $(shell find $(CUDALIB_INC_DIR) -mindepth 0 -type d)
CXX_INC_DIRS 	= $(shell find $(SOLVERS_INC_DIR) -mindepth 1 -type d)
CXX_INC_DIRS 	+= include/ graph/include/ make_trace/include/  ${LIB_INC_DIR} ${CUDA_INC_DIRS} ${SOLVERS_INC_DIR} build/resources/

CXXFLAGS	+=	$(addprefix -I ,$(CXX_INC_DIRS))
NVCCFLAGS	+=	$(addprefix -I ,$(CXX_INC_DIRS))

# Рекурсивный поиск исходных файлов
CXX_SRCS 	:= $(shell find $(CXX_SRC_DIRS) -type f -name '*.cpp')
CUDA_SRCS	:= $(shell find $(CUDA_SRC_DIRS) -type f -name '*.cu')
EXEC_SRCS 	:= $(shell find $(EXE_DIR)/ -type f -name '*.cpp')
TEST_SRCS 	:= $(shell find $(TESTS_DIR)/ -type f -name '*.cpp' ! -name 'off_*.cpp')

# Генерация путей для объектных файлов
CXX_OBJS   	:= $(patsubst %, $(BUILD_OBJ)/%.o, $(CXX_SRCS))
CUDA_OBJS 	:= $(patsubst %, $(BUILD_OBJ)/%.o, $(CUDA_SRCS))
EXECUTABLES := $(patsubst $(EXE_DIR)/%.cpp, $(BUILD_BIN)/%, $(EXEC_SRCS))
TESTS 		:= $(patsubst $(TESTS_DIR)/%.cpp, $(TESTS_BIN)/%, $(TEST_SRCS))

.PHONY: all debug release clean clean-all test

# Запрет автоматического удаления объектных файлов
.PRECIOUS: $(BUILD_OBJ)/%.cpp.o $(BUILD_OBJ)/%.cu.o
.SECONDARY: $(CXX_OBJS) $(CUDA_OBJS)

all: $(EXECUTABLES)
	mkdir -p $(BUILD_DIR)/graph
	mkdir -p $(BUILD_DIR)/illum_geo
	mkdir -p $(BUILD_DIR)/Solve
	mkdir -p $(BUILD_DIR)/trace
	mkdir -p $(BUILD_DIR)/add_dir	

test: $(TESTS)
	mkdir -p $(BUILD_DIR)/graph
	mkdir -p $(BUILD_DIR)/illum_geo
	mkdir -p $(BUILD_DIR)/Solve	

debug:
	@$(MAKE) CONFIG=debug

release:
	@$(MAKE) CONFIG=release

prebuild:	
	$(MAKE) -f $(RESOURCES_DIR)/makefile gen_header_value SRC_FILE=$(ROOT_DIR)/lib/global_consts.h


# Линковка исполняемых файлов тестов
$(TESTS_BIN)/%: $(BUILD_OBJ)/$(TESTS_DIR)/%.cpp.o $(CXX_OBJS) $(CUDA_OBJS)
	@echo "Linking $@"
	@mkdir -p $(@D)
	$(NVCC) -dlink $(addprefix -I ,$(CXX_INC_DIRS)) $(CUDA_OBJS) -o $(notdir $@)_link.o
	$(CXX) -o $@ $^ $(notdir $@)_link.o $(LDFLAGS)
	rm $(notdir $@)_link.o

# Линковка исполняемых файлов
$(BUILD_BIN)/%: $(BUILD_OBJ)/$(EXE_DIR)/%.cpp.o $(CXX_OBJS) $(CUDA_OBJS)
	@echo "Linking $@"
	@mkdir -p $(@D)
	$(NVCC) -dlink $(addprefix -I ,$(CXX_INC_DIRS)) $(CUDA_OBJS) -o $(notdir $@)_link.o
	$(CXX) -o $@ $^ $(notdir $@)_link.o $(LDFLAGS)
	rm $(notdir $@)_link.o

# Компиляция C-файлов
$(BUILD_OBJ)/%.cpp.o: %.cpp
#	@echo "Compiling $<"
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Компиляция CUDA-файлов
$(BUILD_OBJ)/%.cu.o: %.cu
#	@echo "Compiling CUDA $<"
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# Удаление только объектных файлов
clean:
	@echo "Cleaning object files for all configurations"
	@rm -rf $(BUILD_BIN)
	@rm -rf $(BUILD_OBJ)

# Полная очистка (включая исполняемые файлы)
clean-all:
	@echo "Full clean"
	@rm -rf $(BUILD_DIR)

# Стереть логи
log_clean:
	@rm $(BUILD_DIR)/File_Logs*.txt