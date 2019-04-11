cmake \
  -D CMAKE_CXX_COMPILER=mpicxx \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_CXX_FLAGS="-O3 -xHost -ipo -DNDEBUG" \
  -D CMAKE_C_FLAGS="-O3 -xHost -ipo -DNDEBUG" \
  ../
