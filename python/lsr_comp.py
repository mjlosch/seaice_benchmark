import numpy as np
import matplotlib.pyplot as plt
from scipy.io.netcdf import netcdf_file
import cmocean.cm as cmo
import matplotlib.colors as colors
from matplotlib import ticker, cm
import os

rundirs = ['/scratch/users/mlosch/SI_benchmark/4000/runfe00',
           '/scratch/users/mlosch/SI_benchmark/2000/runlsr00',
           '/scratch/users/mlosch/SI_benchmark/2000/runlsr01',
           '/scratch/users/mlosch/SI_benchmark/2000/runfe00']
objects = ('16km', '8km', '4km', '2km')
objects = ('4km,jfnk', '2km, lsr0', '2km, lsr1', '2km, jfnk')
createFig=False

#rundirs = ['./']; createFig=False

t = -1

ieee='b'
accuracy='float64'

def sq(a):
    import numpy as np
    a = np.squeeze(a)
    masked_array=np.ma.masked_where(a==0., a)
    return masked_array

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

functionals = np.zeros((9,len(rundirs)))
fig, ax =plt.subplots(4,3,sharex=True,sharey=True)
fig.set_size_inches(8,11)

for kr, mydir in enumerate(rundirs):
    nc = netcdf_file(os.path.join(mydir,'snapshot.0000000000.t001.nc'),'r')
    h = np.copy(nc.variables['SIheff'][:,:,:,:])
    a = np.copy(nc.variables['SIarea'][:,:,:,:])
    u = np.copy(nc.variables['SIuice'][:,:,:,:])
    v = np.copy(nc.variables['SIvice'][:,:,:,:])
    sh = np.copy(nc.variables['SIshear'][:,:,:,:])
    s1= np.copy(nc.variables['SIsig1'][:])
    s2= np.copy(nc.variables['SIsig2'][:])
    nc.close()

    ncg = netcdf_file(os.path.join(mydir,'grid.t001.nc'),'r')
    xc  = np.copy(ncg.variables['XC'][:,:])
    yc  = np.copy(ncg.variables['YC'][:,:])
    xg  = np.copy(ncg.variables['XG'][:-1,:-1])
    yg  = np.copy(ncg.variables['YG'][:-1,:-1])
    dxg = np.copy(ncg.variables['dxG'][:-1,:])
    dyg = np.copy(ncg.variables['dyG'][:,:-1])
    ncg.close()

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

    functionals[0,kr] = ((totdef*myarea).sum()*ra)
    functionals[1,kr] = ((np.abs(divergnc)*myarea).sum()*ra)
    functionals[2,kr] = ((sheardef*myarea).sum()*ra)
    functionals[3,kr] = ((ke*myarea).sum()*ra)
    functionals[4,kr] = ((hloc*myarea).sum()*ra)
    functionals[5,kr] = u[t,0,:,:-1].min()
    functionals[6,kr] = u[t,0,:,:-1].max()
    functionals[7,kr] = v[t,0,:,:-1].min()
    functionals[8,kr] = v[t,0,:-1,:].max()

    print('area integrals of %s'% mydir)
    print('total deformation     = %e m^2/s'% functionals[0,kr] )
    print('divergent deformation = %e m^2/s'% functionals[1,kr] )
    print('shear deformation     = %e m^2/s'% functionals[2,kr] )
    print('total stress          = %e Nm^2/m'% ((totstress*myarea).sum()*ra))
    print('total kinetic energy  = %e kg m^2/m/s^2'% functionals[3,kr] )
    print('volume of subdomain   = %e m^3'% functionals[4,kr] )

    x = xg*1e-3
    y = yg*1e-3
    csf0=ax[kr,0].pcolormesh(x,y,sq(a[t,0,:,:]),vmin=0.75,vmax=1.,
                             cmap=cmo.ice)
    csf1=ax[kr,1].pcolormesh(x,y,sq(uc[t,0,:,:]),vmin=-.15,vmax=0.15,
                             cmap=cmo.delta)
    csf2=ax[kr,2].pcolormesh(x,y,np.log10(sq(sh[t,0,:,:])),vmin=-8,vmax=-4,
                             cmap=cmo.dense)
    csf = [csf0,csf1,csf2]
    ax[kr,0].set_ylabel(objects[kr])

ax[0,0].set_title('SIarea')
ax[0,1].set_title('SIuice (m/s)')
ax[0,2].set_title('log10(SIshear) (1/s)')
for k, bx in enumerate(ax[-1,:]):
    bx.set_xlabel('(km)')
    pos = bx.get_position()
    cbax = fig.add_axes([pos.x0,pos.y0-0.06,pos.width,0.01])
    if k==0:
        plt.colorbar(csf[k],cax=cbax,orientation='horizontal',extend='min')
    else:
        plt.colorbar(csf[k],cax=cbax,orientation='horizontal',extend='both')

#fig.savefig('lsrcomp')
fig.show()


mytitle = ('total deformation (m$^2$/s)',
           'divergent deformation (m$^2$/s)',
           'shear deformation (m$^2$/s)',
           'kinetic energy (kg m$^2$/s/m$^2$)',
           'volume of subdomain (m$^3$)',
           'min. of x-velocity (m/s)',
           'max. of x-velocity (m/s)',
           'min. of y-velocity (m/s)',
           'max. of y-velocity (m/s)')

if createFig:
    fig, ax =plt.subplots(5,1,sharex=True)
    fig.set_size_inches(6,11)
    y_pos = np.arange(len(objects))
    for k in np.arange(len(ax)):
        ax[k].plot(y_pos, functionals[k,:],'kx-')
        ax[k].set_title(mytitle[k])
        ax[k].grid()

    plt.xticks(y_pos, objects)

#    fig.savefig('functionals')
    plt.show()

    # print values to screen:
    print("%35s , %13s , %13s , %13s , %13s"%
          ("resolution", "16km", "8km", "4km", "2km") )
    for k in np.arange(len(mytitle)):
        print("%35s , %13.6e , %13.6e , %13.6e , %13.6e"%
              (mytitle[k], functionals[k,0], functionals[k,1],
               functionals[k,2], functionals[k,3]))
