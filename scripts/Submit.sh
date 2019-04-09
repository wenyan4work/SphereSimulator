#!/bin/bash

foldername=${PWD##*/} # some bash magic
rm ./*.subscript
cp ./jobsub.slurm ./$foldername.subscript
sbatch -p ccb --qos=ccb  ./$foldername.subscript 
#sbatch  ./$foldername.subscript 
