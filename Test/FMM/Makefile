# mac
# change this to match your system
SFTPATH=/Users/wyan/local

# linux
# SFTPATH=/home_local/software
# EIGEN= $(SFTPATH)/Eigen/Eigen

include $(PVFMM_DIR)/MakeVariables

CXX= mpicxx # must use gcc to link libstfmm3d.a compiled with gfortran
LINK = $(CXX)

CXXFLAGS= $(CXXFLAGS_PVFMM) 
LINKFLAGS= ./libstfmm3d.a -lgfortran $(LDLIBS_PVFMM) $(CXXFLAGS)

# set of libraries to correctly link to the target
USERLIB = -L/Users/wyan/local/lib/ -lmpi
USERINCLUDE = -I./ -I/Users/wyan/local/include/

INCLUDE_DIRS = $(USERINCLUDE)
LIBRARY_DIRS =  
LIBRARIES = $(USERLIB)

# System-specific settings
SHELL = /bin/bash
SYSLIB =	
SIZE =	size

# Files

SRC =   main.cpp fmm.cpp
INC =   fmm.h cmdparser.hpp CustomKernel.hpp
         
# Definitions

EXE =   STFMM3D_INTF.X
OBJ =   $(SRC:.cpp=.o)

# Link rule
$(EXE):	$(OBJ)
	$(LINK) $(OBJ)  -o $(EXE) $(LINKFLAGS) $(SYSLIB) $(LIBRARY_DIRS) $(LIBRARIES)
	$(SIZE) $(EXE)

# Compilation rules
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -c $*.cpp -o $*.o

# Individual dependencies
$(OBJ): $(INC)

all: $(EXE)

clean: 
	rm ./$(OBJ)
	rm ./$(EXE)
