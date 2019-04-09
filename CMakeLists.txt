cmake_minimum_required(VERSION 3.10) # same requirement as trilinos
project(SphereSimulator)

# global settings
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(SFTPATH $ENV{HOME}/local)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)  
# this does not change the compiler
# to use mpi compiler, do cmake -D CMAKE_CXX_COMPILER=mpicxx
# linking executable with mpi
# target_link_libraries(exe PRIVATE OpenMP::OpenMP_CXX)
# target_link_libraries(exe PRIVATE MPI::MPI_CXX)

# external libraries, static link
find_package(TRNG REQUIRED)
find_package(PVFMM REQUIRED)

# external libraries, header only
find_package(Eigen3 REQUIRED)
find_package(MSGPACK REQUIRED)

# internal submodules
set(SIMTOOLBOX ${CMAKE_SOURCE_DIR}/SimToolbox)
set(SCTL ${CMAKE_SOURCE_DIR}/SCTL)
set(STKFMM ${CMAKE_SOURCE_DIR}/STKFMM)

# Get Trilinos as one entity
find_package(Trilinos REQUIRED)
# clean up the too-long list
LIST(REMOVE_DUPLICATES Trilinos_INCLUDE_DIRS)
LIST(REMOVE_DUPLICATES Trilinos_TPL_INCLUDE_DIRS)
LIST(REVERSE Trilinos_LIBRARIES)
LIST(REMOVE_DUPLICATES Trilinos_LIBRARIES)
LIST(REVERSE Trilinos_LIBRARIES)
LIST(REVERSE Trilinos_TPL_LIBRARIES)
LIST(REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)
LIST(REVERSE Trilinos_TPL_LIBRARIES)
# Echo trilinos build info just for fun
MESSAGE("\nFound Trilinos!  Here are the details: ")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("End of Trilinos details\n")

message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
message(STATUS "SFTPATH: ${SFTPATH}")

enable_testing()

set(SRC ${CMAKE_SOURCE_DIR}/SRC)
set(SPHERE ${CMAKE_SOURCE_DIR}/Sphere)
set(MPI ${CMAKE_SOURCE_DIR}/MPI)
set(BIE ${CMAKE_SOURCE_DIR}/BIE)

add_executable(SphereSimulator.X 
${SRC}/Config.cpp ${SRC}/main.cpp ${SRC}/SphereSystem.cpp 
${SPHERE}/Sphere.cpp ${SPHERE}/Shexp.cpp ${BIE}/SphereSTKMobMat.cpp ${BIE}/SphereSTKSHOperator.cpp
${SIMTOOLBOX}/Collision/CollisionSolver.cpp
${SIMTOOLBOX}/Collision/CPSolver.cpp 
${SIMTOOLBOX}/Trilinos/TpetraUtil.cpp
${SIMTOOLBOX}/Util/Base64.cpp
${STKFMM}/STKFMM.cpp
)

target_compile_options(SphereSimulator.X PRIVATE ${OpenMP_CXX_FLAGS})

target_compile_definitions(SphereSimulator.X PRIVATE 
# PVFMM Config
PVFMM_FFTW3_MKL 
# SCTL Configs
# SCTL_MEMDEBUG # Enable memory checks
# SCTL_PROFILE=5 -DSCTL_VERBOSE # Enable profiling
SCTL_HAVE_MPI
SCTL_GLOBAL_MEM_BUFF=512 # Global memory buffer size in MB
SCTL_HAVE_BLAS # use BLAS
SCTL_HAVE_LAPACK # use LAPACK
SCTL_FFTW3_MKL
SCTL_HAVE_FFTW
SCTL_HAVE_FFTWF
)


target_include_directories(SphereSimulator.X PRIVATE
${CMAKE_SOURCE_DIR}
${SRC}
${SCTL}/include
${SIMTOOLBOX}
${STKFMM}
${PVFMM_INCLUDE_DIR}
${TRNG_INCLUDE_DIR}
${EIGEN3_INCLUDE_DIR} 
${MSGPACK_INCLUDE_DIR}
${Trilinos_INCLUDE_DIRS}
$ENV{MKLROOT}/include
$ENV{MKLROOT}/include/fftw
)

target_link_libraries(SphereSimulator.X PRIVATE 
${TRNG_LIBRARY} 
${PVFMM_LIBRARY} 
${Trilinos_LIBRARIES}
${Trilinos_TPL_LIBRARIES}
OpenMP::OpenMP_CXX
MPI::MPI_CXX
)

add_subdirectory(MPI)
add_subdirectory(Sphere)
