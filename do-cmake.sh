cmake \
  -D CMAKE_CXX_COMPILER=mpicxx \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_CXX_FLAGS="-O3 -xHost -ipo -DNDEBUG" \
  -D CMAKE_C_FLAGS="-O3 -xHost -ipo -DNDEBUG" \
  -D SFTPATH="$HOME_MNT/local" \
  -D Eigen3_DIR="$HOME_MNT/local/share/eigen3/cmake/" \
  -D Trilinos_DIR="$HOME_MNT/local/lib/cmake/Trilinos" \
  ../

# gcc flags:
#  -D CMAKE_CXX_FLAGS="-O3 -march=native -DNDEBUG" \
#  -D CMAKE_C_FLAGS="-O3 -march=native -DNDEBUG" \
