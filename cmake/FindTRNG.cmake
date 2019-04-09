# Defines the following variables:
#   - TRNG_LIBRARY
#   - TRNG_INCLUDE_DIR

# Find the header files
find_path(TRNG_INCLUDE_DIR trng/ 
  HINTS ${SFTPATH} $ENV{TRNG_DIR} /usr/local
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with TRNG header.")
MESSAGE("   TRNG_INCLUDE_DIR = ${TRNG_INCLUDE_DIR}")

# Find the library
find_library(TRNG_LIBRARY libtrng4.a 
  HINTS ${SFTPATH} $ENV{TRNG_DIR} /usr/local/
  PATH_SUFFIXES lib/
  NO_DEFAULT_PATH
  DOC "The TRNG library.")
MESSAGE("   TRNG_LIBRARY = ${TRNG_LIBRARY}")
