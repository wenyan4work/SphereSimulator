#!/bin/bash

#NOTE: if you request email, you will get one for each task!
# The stdout (-o) and stderr (-e) files will be written to the submission directory with names that follow the given templates.
#SBATCH -o outrunLog.%j.%A.%a.%N.out 
#SBATCH -e errrunLog.%j.%A.%a.%N.err

#//SBATCH --constraint=broadwell
#SBATCH --constraint=skylake

#SBATCH --nodes=2
#SBATCH --ntasks-per-socket=4
#SBATCH --cpus-per-task=5

#SBATCH --exclusive
#SBATCH --export=ALL

find ./ -type f -name 'outrunLog*' -a ! -name "outrunLog.$SLURM_JOBID.*" -delete 
find ./ -type f -name 'errrunLog*' -a ! -name "errrunLog.$SLURM_JOBID.*" -delete 

date
hostname
echo "Greetings from $SLURM_JOB_NUM_NODES!"
env | grep -i slurm

echo "---------OMP SETTINGS---------"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_DISPLAY_ENV=true
echo $OMP_NUM_THREADS
env | grep -i OMP_
echo "------------------------------"

echo "---------KMP SETTINGS---------"
env | grep -i KMP_
echo "------------------------------"

echo "---------MKL SETTINGS---------"
env | grep -i MKL
echo "------------------------------"

rm -rf ./result/result*-*
rm -rf ./result/ResultAll
find . -type f -name '*.o' -delete
find . -type f -name '*.d' -delete

which mpirun
which mpiexec
which srun
date

echo "hosts:"
echo $SLURM_JOB_NODELIST

export OMP_PROC_BIND=spread
export OMP_PLACES=cores

mpirun --map-by ppr:$SLURM_NTASKS_PER_SOCKET:socket:PE=$SLURM_CPUS_PER_TASK \
       --bind-to hwthread --display-map ./SphereSimulator.X

date
