export OMP_MAX_ACTIVE_LEVELS=1
export OMP_NUM_THREADS=6

mpirun --map-by ppr:1:socket:PE=$OMP_NUM_THREADS \
       --bind-to core --display-map ./SphereSimulator.X > ./outrunLog.out

