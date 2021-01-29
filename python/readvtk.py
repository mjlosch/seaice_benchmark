#
import numpy as np
import matplotlib.pyplot as plt
import meshio
import cmocean.cm as cmo
from myutils import *
from scipy.interpolate import griddata

# Carolin has scaled all her data with these scales
Lscale = 1e6
Tscale = 1e3

filename = './v.00001.vtk'
mesh = meshio.read(filename)

x = mesh.points[:,0]*Lscale
y = mesh.points[:,1]*Lscale
cells = mesh.cells['quad']
vec = mesh.point_data['V'][:,:-1]*Lscale/Tscale

# average to cell centers
x0 = x[cells].mean(axis=1)
y0 = y[cells].mean(axis=1)
u0 = vec[cells,0].mean(axis=1)
v0 = vec[cells,1].mean(axis=1)

# generate a useful grid (similar to MITgcm)
xc,yc=np.meshgrid(np.unique(x0),np.unique(y0))
# and interpolate (map) to it
uc = griddata((x0,y0), u0, (xc,yc), method='nearest')
vc = griddata((x0,y0), v0, (xc,yc), method='nearest')

# generate the corner grid
xg,yg=np.meshgrid(np.unique(x),np.unique(y))
# and interpolate (map) to it
ug = griddata((x,y), vec[:,0], (xg,yg), method='nearest')
vg = griddata((x,y), vec[:,1], (xg,yg), method='nearest')

# some testplot
fig=plt.figure()
plt.pcolormesh(xg,yg,sq(uc))
plt.colorbar()
plt.contour(xg,yg,sq(ug),20,colors='k')

plt.show()
