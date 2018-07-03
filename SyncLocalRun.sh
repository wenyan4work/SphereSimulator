#! /bin/bash
RUNPATH=/local1/RunLocal/SphereVWX

# copy to cluster run path
rsync -avz --delete-after --exclude='*.d' --exclude='*.o' --exclude='*.png' --exclude='*.vtp' --exclude='*.xyz' --exclude='*.vtk' --exclude=".*" ./* $RUNPATH/
