# CCB014
SFTPATH=/mnt/home/wyan/local

# Mac 
#SFTPATH=/Users/wyan/local

include $(SFTPATH)/include/Makefile.export.Trilinos
include $(PVFMM_DIR)/MakeVariables
TRNG = $(SFTPATH)/include/trng
EIGEN= $(SFTPATH)/include/eigen3
PVFMM= $(SFTPATH)/include/pvfmm
FDPS=./FDPS_Header


USERLIB = -lm -ltrng4
USERINCLUDE =  -I$(FDPS) -I$(TRNG)/include -I$(EIGEN) -I$(MKLROOT)/include/fftw -I$(PVFMM)

CXX=mpicxx

LINK = $(CXX)

# optimized
CXXFLAGS= $(Trilinos_CXX_COMPILER_FLAGS) $(CXXFLAGS_PVFMM) -Wa,-q #-DDEBUGLCPCOL


LINKFLAGS= $(CXXFLAGS) $(Trilinos_EXTRA_LD_FLAGS)  $(LDLIBS_PVFMM) -lm -ldl

# debug
#CXXFLAGS = -std=c++11 -O0 -fopenmp=libiomp5 -g #-DDEBUGLCPCOL #-DZDDDEBUG -DIFPACKDEBUG #-DMYDEBUGINFO -DPARTICLE_SIMULATOR_DEBUG_PRINT -DHYRDRODEBUG
#LINKFLAGS= $(CXXFLAGS) $(Trilinos_EXTRA_LD_FLAGS)  $(LDLIBS_PVFMM) -lm -ldl

# almost always yes
WITHOPENMP ?= yes
WITHMPI ?= yes

ifeq ($(WITHOPENMP), yes)
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
endif

ifeq ($(WITHMPI), yes)
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
endif
