import numpy as np
from myutils import *
#from matplotlib.mlab import find
import matplotlib.pyplot as plt
from scipy.io.netcdf import netcdf_file
import cmocean.cm as cmo
import matplotlib.colors as colors
from matplotlib import ticker, cm

ieee='b'
accuracy='float64'

nc = netcdf_file('../test/snapshot.0000000000.t001.nc','r')
h = np.copy(nc.variables['SIheff'][:,:,:,:])
a = np.copy(nc.variables['SIarea'][:,:,:,:])
u = np.copy(nc.variables['SIuice'][:,:,:,:])
v = np.copy(nc.variables['SIvice'][:,:,:,:])
s1= np.copy(nc.variables['SIsig1'][:])
s2= np.copy(nc.variables['SIsig2'][:])
nc.close()

ncg = netcdf_file('../test/grid.t001.nc','r')
xg  = np.copy(ncg.variables['XG'][:-1,:-1])
yg  = np.copy(ncg.variables['YG'][:-1,:-1])
dxg = np.copy(ncg.variables['dxG'][:-1,:])
dyg = np.copy(ncg.variables['dyG'][:,:-1])
ncg.close()

t = 1

def strainrates( ui, vi, dxg, dyg ):
    dxv = 0.5*(dxg+np.roll(dxg,1,1))
    dyu = 0.5*(dyg+np.roll(dyg,1,0))
    dxf = 0.5*(dxg+np.roll(dxg,-1,0))
    dyf = 0.5*(dyg+np.roll(dyg,-1,1))
    e11=(np.roll(ui,-1,1)-ui)/dxf
    e22=(np.roll(vi,-1,0)-vi)/dyf
    e12sql=0.25*( (ui-np.roll(ui,1,0))/dyu + (vi-np.roll(vi,1,1))/dxv )**2
    e12sq = 0.25*(e12sql + np.roll(e12sql,-1,0) + np.roll(e12sql,-1,1)
                  + np.roll(np.roll(e12sql,-1,0),-1,1))
    return e11, e22, e12sq

def deformation( e1, e2, e12sq ):

    divergnc = e1
    sheardef = np.sqrt(e2**2 + 4.*e12sq)

    return divergnc, sheardef

def kinEnergy(h,ui,vi):

    KE = 900 * 0.5 *h * ( 0.5*(ui[:,:-1]**2 + ui[:,1:]**2)
                        + 0.5*(vi[:-1,:]**2 + vi[1:,:]**2) )
    return KE

uc = 0.5*(u[:,:,:,:-1]+u[:,:,:,1:])
vc = 0.5*(v[:,:,:-1,:]+v[:,:,1:,:])

e11,e22,e12sq = strainrates(u[t,0,:,:-1],v[t,0,:-1,:],dxg,dyg)

divergnc, sheardef = deformation(e11+e22,e11-e22,e12sq)

totdef = np.sqrt(divergnc**2+sheardef**2)

ke = kinEnergy(h[t,0,:,:],u[t,0,:,:],v[t,0,:,:])

totstress = np.sqrt((s1+s2)**2 + (s1-s2)**2)[t,0,:,:]

myarea = dxg*dyg
myarea[:,-1] = 0.
myarea[-1,:] = 0.
ra = 1./(myarea.sum())

ra = 1.
#myarea = 1.

hloc = np.copy(h[t,0,:,:])
hloc[np.logical_or(xg<375e3,yg<375e3)] = 0

print('area integrals of')
print('total deformation     = %e m^2/s'% ((totdef*myarea).sum()*ra) )
print('divergent deformation = %e m^2/s'% ((np.abs(divergnc)*myarea).sum()*ra) )
print('shear deformation     = %e m^2/s'% ((sheardef*myarea).sum()*ra) )
print('total stress          = %e Nm^2/m'% ((totstress*myarea).sum()*ra))
print('total kinetic energy  = %e kg m^2/m/s^2'% ((ke*myarea).sum()*ra))
print('volume of subdomain   = %e m^3'% ((hloc*myarea).sum()*ra))

plt.clf()
plt.pcolormesh(xg/1000,yg/1000,sq(totdef[:-1,:-1]),
               norm=colors.LogNorm(vmin=sq(totdef).min(), vmax=totdef.max()),
               cmap=cmo.balance)
plt.colorbar()

plt.show()
