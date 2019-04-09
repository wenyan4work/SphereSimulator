# Example of how to compile and link
# target_compile_options
# target_include_directories
# target_link_libraries

add_executable(Buffer_test Buffer_test.cpp ${SPHERE}/Sphere.cpp ${SPHERE}/Shexp.cpp ${SIMTOOLBOX}/Util/Base64.cpp)
target_compile_options(Buffer_test PRIVATE ${OpenMP_CXX_FLAGS})
target_compile_definitions(Buffer_test PRIVATE 
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

target_include_directories(Buffer_test PRIVATE
${CMAKE_SOURCE_DIR}
${SRC}
${SCTL}/include
${SIMTOOLBOX}
${EIGEN3_INCLUDE_DIR} 
${MSGPACK_INCLUDE_DIR}
$ENV{MKLROOT}/include
$ENV{MKLROOT}/include/fftw
)

target_link_libraries(Buffer_test PRIVATE 
${Trilinos_LIBRARIES}
${Trilinos_TPL_LIBRARIES}
OpenMP::OpenMP_CXX
MPI::MPI_CXX
)

add_executable(InteractionManager_test InteractionManager_test.cpp)
target_compile_options(InteractionManager_test PRIVATE ${OpenMP_CXX_FLAGS})
target_compile_definitions(InteractionManager_test PRIVATE 
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

target_include_directories(InteractionManager_test PRIVATE
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

target_link_libraries(InteractionManager_test PRIVATE 
${TRNG_LIBRARY} 
${PVFMM_LIBRARY} 
${Trilinos_LIBRARIES}
${Trilinos_TPL_LIBRARIES}
OpenMP::OpenMP_CXX
MPI::MPI_CXX
)


add_executable(SyncManager_test SyncManager_test.cpp)
target_compile_options(SyncManager_test PRIVATE ${OpenMP_CXX_FLAGS})
target_include_directories(SyncManager_test PRIVATE ${MSGPACK_INCLUDE_DIR})
target_link_libraries(SyncManager_test PRIVATE OpenMP::OpenMP_CXX MPI::MPI_CXX)