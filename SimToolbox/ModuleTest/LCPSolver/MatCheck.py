import numpy as np
import scipy as sp
import scipy.io as so
import matplotlib as mpl
import matplotlib.pyplot as plt

import os

os.system('export OMP_NUM_THREADS=2 && mpiexec -n 6 ./TestLCPSolver.X 200 > ./test.log ')

# show Amat
Amat = so.mmread('Amat_TCMAT.mtx').todense()
bvec = so.mmread('bvec_TV.mtx')

fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(Amat, interpolation='nearest', norm=mpl.colors.LogNorm())
fig.colorbar(cax)
plt.savefig('Amat.png', dpi=300)
