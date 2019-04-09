SFTPATH:=$(HOME)/local

# inherit flags, includes and libs from Trilinos and pvfmm
include $(SFTPATH)/include/Makefile.export.Trilinos
include $(PVFMM_DIR)/MakeVariables

# internal includes
SCTL := ${CURDIR}/SCTL/include
SIMTOOLBOX := ${CURDIR}/SimToolbox
STKFMM := ${CURDIR}/STKFMM

# external libraries
TRNG  := $(SFTPATH)/include/trng
EIGEN := $(SFTPATH)/include/eigen3
PVFMM := $(SFTPATH)/include/pvfmm

CXX= mpicxx
LINK= $(CXX)

# optimized
CXXFLAGS= $(CXXFLAGS_PVFMM) -ipo -xHost
LINKFLAGS= $(CXXFLAGS) $(Trilinos_EXTRA_LD_FLAGS)

# remove some flags for debugging

# debug
DEBUGMODE:= no

# debug flags
# CXXFLAGS += -DFMMTIMING 
# CXXFLAGS += -DFMMDEBUG
# CXXFLAGS += -DDEBUGLCPCOL 
# CXXFLAGS += -DZDDDEBUG 
# CXXFLAGS += -DIFPACKDEBUG 
# CXXFLAGS += -DMYDEBUGINFO 
# CXXFLAGS += -DHYRDRODEBUG

ifeq ($(DEBUGMODE), yes)
CXXFLAGS:= $(subst -O3, ,$(CXXFLAGS))
LINKFLAGS:= $(subst -O3, ,$(LINKFLAGS))
CXXFLAGS := $(CXXFLAGS) -O0 -g
LINKFLAGS := $(LINKFLAGS) -O0 -g
else
CXXFLAGS:= $(CXXFLAGS) -DNDEBUG
LINKFLAGS:= $(LINKFLAGS) -DNDEBUG
endif

# SCTL Configs

# almost always yes
WITHMPI ?= yes

ifeq ($(WITHMPI), yes)
CXXFLAGS += -DSCTL_HAVE_MPI
endif

# CXXFLAGS += -DSCTL_MEMDEBUG # Enable memory checks
CXXFLAGS += -DSCTL_GLOBAL_MEM_BUFF=512 # Global memory buffer size in MB
# CXXFLAGS += -DSCTL_PROFILE=5 -DSCTL_VERBOSE # Enable profiling
# CXXFLAGS += -DSCTL_QUAD_T=__float128 # Enable quadruple precision

CXXFLAGS += -DSCTL_HAVE_BLAS # use BLAS
CXXFLAGS += -DSCTL_HAVE_LAPACK # use LAPACK
CXXFLAGS += -DSCTL_HAVE_FFTW
CXXFLAGS += -DSCTL_HAVE_FFTWF

# CXXFLAGS += -mkl -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK # use MKL BLAS and LAPACK
# CXXFLAGS += -lblas -DSCTL_HAVE_BLAS # use BLAS
# CXXFLAGS += -llapack -DSCTL_HAVE_LAPACK # use LAPACK
# CXXFLAGS += -lfftw3 -DSCTL_HAVE_FFTW
# CXXFLAGS += -lfftw3f -DSCTL_HAVE_FFTWF
# CXXFLAGS += -lfftw3l -DSCTL_HAVE_FFTWL
