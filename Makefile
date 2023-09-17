# defining working directories
SRCDIR          = src src/fun_src
INCLUDESDIR     = include include/fun_inc
BUILDDIR        = build
OBJDIR          = $(BUILDDIR)/objs

# specify the list of directories in which the search should be performed.
vpath %.cpp 	$(SRCDIR) 
vpath %.h 		$(INCLUDESDIR)
vpath %.o 		$(OBJDIR)
vpath %.d 		$(BUILDDIR)

.SUFFIXES:					# Delete the default suffixes
.SUFFIXES: .cpp .h .o .d	# Define our suffix list


SRCS 			= $(foreach dir,$(SRCDIR),$(wildcard $(dir)/*.cpp))
OBJS            = $(patsubst %.cpp, %.o, $(notdir $(SRCS)))
DEP_FILES       = $(patsubst %.o, %.d, $(OBJS))
INCLUDE_DIRS    = $(addprefix -I ,$(INCLUDESDIR))

# configuring the compiler
CXX             = mpic++
CPPFLAGS        = 
CXXFLAGS        = -g #-Wall -Wextra -std=c++11
PROGRAM         = run
#TARGET_ARCH		=

COMPILE.cpp     = $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDE_DIRS) $(TARGET_ARCH)

all: $(PROGRAM)

$(PROGRAM): %: $(OBJS)
	$(LINK.cpp) $(INCLUDE_DIRS) $(addprefix $(OBJDIR)/, $^) $(LOADLIBES) $(LDLIBS) -o $@


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