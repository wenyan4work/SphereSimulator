# CCB014
SFTPATH=/mnt/home/wyan/local

# Mac 
#SFTPATH=/Users/wyan/local

include $(SFTPATH)/include/Makefile.export.Trilinos
include $(PVFMM_DIR)/MakeVariables
TRNG = $(SFTPATH)/include/trng
EIGEN= $(SFTPATH)/include/eigen3
PVFMM= $(SFTPATH)/include/pvfmm

USERINCLUDE = -I$(TRNG)/include -I$(EIGEN) -I$(PVFMM) # -I$(MKLROOT)/include/fftw 
USERLIB_DIRS = -L$(SFTPATH)/lib
USERLIBS = -ltrng4

CXX=mpicxx

LINK = $(CXX)

# optimized
CXXFLAGS= $(CXXFLAGS_PVFMM) #-DDEBUGLCPCOL
LINKFLAGS= $(CXXFLAGS) $(LDLIBS_PVFMM) $(Trilinos_EXTRA_LD_FLAGS) #-lm -ldl

#CXXFLAGS= $(Trilinos_CXX_COMPILER_FLAGS) $(CXXFLAGS_PVFMM) #-DDEBUGLCPCOL
#LINKFLAGS= $(CXXFLAGS) $(Trilinos_EXTRA_LD_FLAGS)  $(LDLIBS_PVFMM) #-lm -ldl

# debug
#CXXFLAGS = -std=c++11 -O0 -fopenmp=libiomp5 -g #-DDEBUGLCPCOL #-DZDDDEBUG -DIFPACKDEBUG #-DMYDEBUGINFO -DPARTICLE_SIMULATOR_DEBUG_PRINT -DHYRDRODEBUG
#LINKFLAGS= $(CXXFLAGS) $(Trilinos_EXTRA_LD_FLAGS)  $(LDLIBS_PVFMM) -lm -ldl
