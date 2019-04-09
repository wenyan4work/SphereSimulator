SFTPATH:=$(HOME)/local

# inherit flags, includes and libs from Trilinos and pvfmm
include $(SFTPATH)/include/Makefile.export.Trilinos
include $(PVFMM_DIR)/MakeVariables

# internal includes
SCTL = ${CURDIR}/SCTL/include
SIMTOOLBOX = ${CURDIR}/SimToolbox
STKFMM = ${CURDIR}/STKFMM

# external libraries
# static link
TRNG_DIR  = $(SFTPATH)/include/trng
TRNG_LIB  = $(SFTPATH)/lib/libtrng4.a
PVFMM_DIR = $(SFTPATH)/include/pvfmm
PVFMM_LIB = $(SFTPATH)/lib/pvfmm/libpvfmm.a

# header only libraries
EIGEN  = $(SFTPATH)/include/eigen3
MSGPACK = $(SFTPATH)/include/

CXX  = mpicxx
LINK = $(CXX)

# optimized
CXXFLAGS  := -O3 -std=c++14 -qopenmp -qno-offload -ipo -xHost
LINKFLAGS := $(CXXFLAGS) $(Trilinos_EXTRA_LD_FLAGS)

# debug
DEBUGMODE = no

# debug flags
# CXXFLAGS += -DFMMTIMING 
# CXXFLAGS += -DFMMDEBUG
# CXXFLAGS += -DDEBUGLCPCOL 
# CXXFLAGS += -DZDDDEBUG 
# CXXFLAGS += -DIFPACKDEBUG 
# CXXFLAGS += -DMYDEBUGINFO 
# CXXFLAGS += -DHYRDRODEBUG

# remove some flags for debugging
ifeq ($(DEBUGMODE), yes)
CXXFLAGS  = $(subst -O3, ,$(CXXFLAGS))
LINKFLAGS = $(subst -O3, ,$(LINKFLAGS))
CXXFLAGS  = $(CXXFLAGS) -O0 -g
LINKFLAGS = $(LINKFLAGS) -O0 -g
else
CXXFLAGS  = $(CXXFLAGS) -DNDEBUG
LINKFLAGS = $(LINKFLAGS) -DNDEBUG
endif

# PVFMM Configs
CXXFLAGS += -DPVFMM_FFTW3_MKL

# SCTL Configs
# CXXFLAGS += -DSCTL_MEMDEBUG # Enable memory checks
# CXXFLAGS += -DSCTL_PROFILE=5 -DSCTL_VERBOSE # Enable profiling
CXXFLAGS += -DSCTL_HAVE_MPI
CXXFLAGS += -DSCTL_GLOBAL_MEM_BUFF=512 # Global memory buffer size in MB
CXXFLAGS += -DSCTL_HAVE_BLAS # use BLAS
CXXFLAGS += -DSCTL_HAVE_LAPACK # use LAPACK
CXXFLAGS += -DSCTL_FFTW3_MKL
CXXFLAGS += -DSCTL_HAVE_FFTW
CXXFLAGS += -DSCTL_HAVE_FFTWF

# set of libraries to correctly link to the target
USERINC_DIRS = -I$(TRNG_DIR) -I$(PVFMM_DIR) -I$(EIGEN) -I$(MSGPACK) \
               -I$(SCTL) -I$(SIMTOOLBOX) -I$(CURDIR) \
			     $(FFTW_INCLUDES_PVFMM)
USERLIB_DIRS = -L$(SFTPATH)/lib
USERLIBS = $(TRNG_LIB) $(PVFMM_LIB)

INCLUDE_DIRS = $(Trilinos_INCLUDE_DIRS) $(Trilinos_TPL_INCLUDE_DIRS) $(USERINC_DIRS)
LIBRARY_DIRS = $(Trilinos_LIBRARY_DIRS) $(Trilinos_TPL_LIBRARY_DIRS) $(USERLIB_DIRS)
LIBRARIES = $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) $(USERLIBS)

# System-specific settings
SHELL = /bin/sh
SYSLIB =	
SIZE =	size
