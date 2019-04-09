# Defines the following variables:
#   - PVFMM_LIBRARY
#   - PVFMM_INCLUDE_DIR

# Find the header files
find_path(PVFMM_INCLUDE_DIR pvfmm.hpp
  HINTS ${SFTPATH} $ENV{PVFMM_DIR} /usr/local
  PATH_SUFFIXES include/pvfmm 
  NO_DEFAULT_PATH
  DOC "Directory with PVFMM header.")
MESSAGE("   PVFMM_INCLUDE_DIR = ${PVFMM_INCLUDE_DIR}")

# Find the library
find_library(PVFMM_LIBRARY pvfmm/libpvfmm.a 
  HINTS ${SFTPATH} $ENV{PVFMM_DIR} /usr/local/
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  DOC "The PVFMM library.")
MESSAGE("   PVFMM_LIBRARY = ${PVFMM_LIBRARY}")