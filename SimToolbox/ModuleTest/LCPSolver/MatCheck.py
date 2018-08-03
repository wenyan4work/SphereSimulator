import numpy as np
import scipy as sp
import scipy.io as so
import matplotlib.pyplot as plt

Amat=so.mmread('Amat_TCMAT.mtx').todense()
bvec=so.mmread('bvec_TV.mtx')

plt.matshow(Amat)
plt.show()
