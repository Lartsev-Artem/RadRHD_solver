#export BUILDDIR
RESDIR  = $(ROOT_DIR)/$(BUILD_DIR)/resources
EXE 	= $(RESDIR)/parser_header
#SRC_FILE - arg make

CUR_DIR :=	$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
#CUR_DIR_NAME := $(dir $(lastword $(CUR_DIR))) 		#путь в формате path/
#CUR_DIR_NAME := $(patsubst %/,%,$(CUR_DIR_NAME))	#путь в формате path
#ROOT_DIR	 := $(patsubst %$(LOOC_DIR),%,$(CUR_DIR))	#путь главного 

all: gen_header_value

gen_header_value: 	
	mkdir -p $(BUILD_DIR)	
	mkdir -p $(RESDIR)
	$(CXX) -std=c++11 -o $(EXE) $(CUR_DIR)/parser_header.cpp -I $(ROOT_DIR)/lib/calculator/ -I $(ROOT_DIR)/lib/output_interface/
	$(EXE) $(SRC_FILE) $(RESDIR)/
	rm $(EXE)
