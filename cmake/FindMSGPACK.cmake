# Defines the following variables:
#   - MSGPACK_INCLUDE_DIR

# Find the header files
find_path(
  MSGPACK_INCLUDE_DIR msgpack.hpp
  HINTS ${SFTPATH} $ENV{PVFMM_DIR} /usr/local
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with PVFMM header.")
message("   MSGPACK_INCLUDE_DIR= ${MSGPACK_INCLUDE_DIR}")
