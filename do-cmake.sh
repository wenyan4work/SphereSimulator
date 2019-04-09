cmake \
  -D CMAKE_CXX_COMPILER=mpicxx \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_CXX_FLAGS="-qno-offload" \
  -D CMAKE_C_FLAGS="-qno-offload" \
  ../
