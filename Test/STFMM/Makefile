# mac
# change this to match your system
# SFTPATH=/Users/wyan/local

# linux
SFTPATH=/home_local/software
EIGEN= $(SFTPATH)/include/eigen3

include $(PVFMM_DIR)/MakeVariables

CXX= mpicxx 
# must use gcc/g++ to link libstfmm3d.a compiled with gfortran
# and use icc/icpc to link libstfmm3d.a compiled with ifort
LINK = mpifc -cxxlib -nofor_main

CXXFLAGS= $(CXXFLAGS_PVFMM) 
LINKFLAGS= $(LDLIBS_PVFMM) $(CXXFLAGS)

# add -lgfortran for gcc linking gfortran
# use ifort to link ifort and icc

# set of libraries to correctly link to the target
# USERLIB = -L/Users/wyan/local/lib/ -lmpi
# USERINCLUDE = -I./ -I/Users/wyan/local/include/

INCLUDE_DIRS = $(USERINCLUDE)
LIBRARY_DIRS =  
LIBRARIES = $(USERLIB)

# System-specific settings
SHELL = /bin/bash
SYSLIB =	
SIZE =	size

# Files

SRC =   main.cpp fmm.cpp
INC =   fmm.h cmdparser.hpp StokesCustomKernel.hpp
         
# Definitions

EXE =   STFMM3D_INTF.X
OBJ :=   $(SRC:.cpp=.o)

FSRC := $(wildcard ./stfmmlibsrc/*.f)
FOBJ := $(FSRC:.f=.o)

# Link rule
$(EXE):	$(OBJ) $(FOBJ)
	$(LINK) $(OBJ) $(FOBJ)  -o $(EXE) $(LINKFLAGS) $(SYSLIB) $(LIBRARY_DIRS) $(LIBRARIES)
	$(SIZE) $(EXE)

# Compilation rules
.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -c $*.cpp -o $*.o

# Individual dependencies
$(OBJ): $(INC)

$(FOBJ):
	cd ./stfmmlibsrc && make

all: $(EXE)

clean: 
	rm ./$(OBJ)
	rm ./stfmmlibsrc/$(FOBJ)
	rm ./$(EXE)
