#! /bin/bash
RUNPATH=/local1/RunLocal/SphereSimulator

# copy to cluster run path
rsync -avz --delete-after --exclude='*.d' --exclude='*.o' \
           --exclude='*.png' --exclude='*.xyz' --exclude='*.mtx' \
           --exclude='*.log' --exclude='*.vtr' --exclude='*.pvtr' \
           --exclude='*.vtk' --exclude='*.vtp' --exclude='*.vtu' \
           --exclude='*.pvtk' --exclude='*.pvtp' --exclude='*.pvtu' \
           --exclude=".*" ./ $RUNPATH