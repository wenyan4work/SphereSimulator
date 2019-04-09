# generate random points in a circle

import numpy as np

nparticle = 20000
radius = 180
height = 0.1

# make a simple unit circle
theta = np.linspace(0, 2*np.pi, nparticle)
a, b = 1 * np.cos(theta), 1 * np.sin(theta)

# generate uniform random points in circle
theta = np.random.rand(nparticle) * (2 * np.pi)
r = np.sqrt(np.random.rand(nparticle))*radius
x, y = r * np.cos(theta), r * np.sin(theta)
z = np.random.rand(nparticle)*height
size = np.random.lognormal(0.0, 0.2, nparticle)

for i in range(len(size)):
    if size[i] < 0.5:
        size[i] = 0.5
    elif size[i] > 2:
        size[i] = 2
    size[i]=1.0

print(nparticle)
print("min radius = "+str(np.min(size)))

for i in range(nparticle):
    print('S '+str(i)+' '+str(size[i])+' '+str(x[i])+' '+str(y[i])+' '+str(z[i]))
