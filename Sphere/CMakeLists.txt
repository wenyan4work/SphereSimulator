# Example of how to compile and link
# target_compile_options
# target_include_directories
# target_link_libraries

add_executable(Sphere_test Sphere_test.cpp Sphere.cpp Shexp.cpp ${SIMTOOLBOX}/Util/Base64.cpp)
target_compile_options(Sphere_test PRIVATE ${OpenMP_CXX_FLAGS})
target_compile_definitions(Sphere_test PRIVATE 
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

target_include_directories(Sphere_test PRIVATE
${CMAKE_SOURCE_DIR}
${SCTL}/include
${SIMTOOLBOX}
${EIGEN3_INCLUDE_DIR} 
${MSGPACK_INCLUDE_DIR}
$ENV{MKLROOT}/include
$ENV{MKLROOT}/include/fftw
)

target_link_libraries(Sphere_test PRIVATE 
${Trilinos_LIBRARIES}
${Trilinos_TPL_LIBRARIES}
OpenMP::OpenMP_CXX
MPI::MPI_CXX
)

add_executable(Shexp_test Shexp_test.cpp Shexp.cpp ${SIMTOOLBOX}/Util/Base64.cpp ${SIMTOOLBOX}/Util/PointDistribution.cpp)
target_compile_options(Shexp_test PRIVATE ${OpenMP_CXX_FLAGS})
target_compile_definitions(Shexp_test PRIVATE 
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

target_include_directories(Shexp_test PRIVATE
${CMAKE_SOURCE_DIR}
${SCTL}/include
${SIMTOOLBOX}
${EIGEN3_INCLUDE_DIR} 
${MSGPACK_INCLUDE_DIR}
$ENV{MKLROOT}/include
$ENV{MKLROOT}/include/fftw
)

target_link_libraries(Shexp_test PRIVATE 
${Trilinos_LIBRARIES}
${Trilinos_TPL_LIBRARIES}
OpenMP::OpenMP_CXX
MPI::MPI_CXX
)