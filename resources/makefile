#export BUILDDIR
RESDIR  = $(BUILDDIR)/resources
EXE 	= $(RESDIR)/run
#SRC_FILE - arg make

ROOT_DIR :=	$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

PREFIX = log

all: gen_header_value

gen_header_value: 
	mkdir -p $(RESDIR)
	$(CXX) -std=c++17 -o $(EXE) $(ROOT_DIR)/generate_header.cpp
	./$(EXE) $(SRC_FILE) $(RESDIR)/$(PREFIX)_$(notdir $(SRC_FILE))
	rm $(EXE)

		