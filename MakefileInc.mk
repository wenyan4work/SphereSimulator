# CCB014
SFTPATH=/mnt/home/wyan/local

# Mac 
#SFTPATH=/Users/wyan/local

include $(SFTPATH)/include/Makefile.export.Trilinos
include $(PVFMM_DIR)/MakeVariables
TRNG = $(SFTPATH)/include/trng
EIGEN= $(SFTPATH)/include/eigen3
PVFMM= $(SFTPATH)/include/pvfmm

USERINCLUDE = -I$(TRNG)/include -I$(EIGEN) # -I$(PVFMM) # -I$(MKLROOT)/include/fftw 
USERLIB_DIRS = -L$(SFTPATH)/lib
USERLIBS = -ltrng4

CXX= mpicxx
LINK= $(CXX)

# optimized
CXXFLAGS= $(CXXFLAGS_PVFMM)
LINKFLAGS= $(CXXFLAGS) $(LDLIBS_PVFMM) $(Trilinos_EXTRA_LD_FLAGS) #-lm -ldl

# remove some flags for debugging
# if Trilinos and pvfmm are compiled with ipo, removing this may cause linking failures

# debug
DEBUGMODE:= yes

# debug flags
# CXXFLAGS += -DFMMTIMING 
# CXXFLAGS += -DDEBUGLCPCOL 
# CXXFLAGS += -DZDDDEBUG 
# CXXFLAGS += -DIFPACKDEBUG 
# CXXFLAGS += -DMYDEBUGINFO 
# CXXFLAGS += -DHYRDRODEBUG
# CXXFLAGS += -DDEBUGLCPCOL

ifeq ($(DEBUGMODE), yes)
CXXFLAGS:= $(subst -O3, ,$(CXXFLAGS))
LINKFLAGS:= $(subst -O3, ,$(LINKFLAGS))
CXXFLAGS := $(CXXFLAGS) -O0 -g
LINKFLAGS := $(LINKFLAGS) -O0 -g
endif

# almost always yes
WITHMPI ?= yes

ifeq ($(WITHMPI), yes)
CXXFLAGS += -DSCTL_HAVE_MPI
endif
