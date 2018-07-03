#! /bin/bash
CLUSTERRUNPATH=/mnt/ceph/users/wyan/RunCluster/SphereVWX

# copy to cluster run path
rsync -avz --delete-after --exclude='*.d' --exclude='*.o' --exclude='*.X' --exclude='*.png' --exclude='*.xyz' --exclude='*.vtp' \
           --exclude='*.vtk' --exclude='*.vtp' --exclude='*.vtu' \
           --exclude='*.pvtk' --exclude='*.pvtp' --exclude='*.pvtu' \
           --exclude=".*" ./* $CLUSTERRUNPATH/
