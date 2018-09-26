#!/bin/bash

#NOTE: if you request email, you will get one for each task!
# The stdout (-o) and stderr (-e) files will be written to the submission directory with names that follow the given templates.
#SBATCH -o outrunLog.%j.%A.%a.%N.out 
#SBATCH -e errrunLog.%j.%A.%a.%N.err

#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --export=ALL

# for one mpi rank per socket
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=14

find ./ -type f -name 'outrunLog*' -a ! -name "outrunLog.$SLURM_JOBID.*" -delete 
find ./ -type f -name 'errrunLog*' -a ! -name "errrunLog.$SLURM_JOBID.*" -delete 

date
hostname
echo "Greetings from $SLURM_JOB_NUM_NODES!"
env | grep -i slurm

echo "---------OMP SETTINGS---------"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo $OMP_NUM_THREADS
env | grep -i OMP
echo "------------------------------"

echo "---------KMP SETTINGS---------"
env | grep -i KMP
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

# manually assign running nodes
# by default each rank binds to one socket, do not specify extra -bind-to or -map-to parameters
# if -bind-to or -map-to parameters exist, job is suspended every 30 seconds for unknown reasons
echo "hosts:"
echo $SLURM_JOB_NODELIST
mpiexec -host $(scontrol show hostname $SLURM_NODELIST | tr '\n' ',') -perhost $SLURM_TASKS_PER_NODE ./SphereSimulator.X


date
