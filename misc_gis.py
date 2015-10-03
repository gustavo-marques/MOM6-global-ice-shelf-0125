# Gustavo Marques
#
#  These are general functions used
#  in the other scripts
#

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from midas.rectgrid import *
import netCDF4
import os

def plt_latlon(lon,lat,var,depth,varname,varunits,levs,date,ind,mymap=plt.cm.jet):
    if not os.path.exists(varname):
        os.system('mkdir '+ varname)

    fig = plt.figure()
    ax = fig.add_subplot(111,axisbg='gray')
    cm=ax.contourf(lon,lat,var,levs,cmap=mymap,extend='both')
    ax.contour(lon,lat,D,5,colors='w',linewidths=2)
    cbar=plt.colorbar(cm,orientation='horizontal',
                     ticks=[levs.min(),(np.abs(levs.max())-np.abs(levs.min()))/2.0,levs.max()])
    cbar.set_label(r'%s [%s]' %(str(varname,varunits)))
    #ax.set_aspect('auto')
    #cs.cmap.set_under('b')
    #cs.cmap.set_over('r')
    #cs.set_clim(levs[0], levs[-1])
    ax.set_title('%s - %s' %(str(varname,date)))
    ax.set_xlabel('Lon, deg E')
    ax.set_ylabel('Lat, deg N')
    plt.grid()
    plt.tight_layout()
    s = str("plt.savefig('%s/%s-%04d.png',bbox_inches='tight')"% (varname,exp,ind))
    eval(s)
    plt.close('all')
    return

def get_grid():
    """
    Get total grid info for the globe
    Usage: grid = get_grid()
    """
    grid_path='/net2/mjh/ipynb/GIS_0125/'
    sgrid=supergrid(file=grid_path+'ocean_hgrid.nc',cyclic_x=True,tripolar_n=True)
    grid=quadmesh(supergrid=sgrid)
    grid.lonh=grid.x_T[grid.jm/2,:]
    grid.lath=grid.y_T[:,grid.im/4]
    grid.lonq=grid.x_T_bounds[grid.jm/2,:]
    grid.latq=grid.y_T_bounds[:,grid.im/4]
    D=netCDF4.Dataset(grid_path+'topog_v4.nc').variables['depth'][:,:]
    grid.D=D
    grid.wet[grid.D<10.]=0.
    return grid

def get_region_grid(lon1,lon2,lat1,lat2):
    """
    Get grid info for a specific region, specified by lat1,lat2,lon1,lon2 
    Usage: grid,region = get_region_grid(lat1,lat2,lon1,lon2)
    """

    grid_path='/net2/mjh/ipynb/GIS_0125/'
    sgrid=supergrid(file=grid_path+'ocean_hgrid.nc',cyclic_x=True,tripolar_n=True)
    grid=quadmesh(supergrid=sgrid)
    grid.lonh=grid.x_T[grid.jm/2,:]
    grid.lath=grid.y_T[:,grid.im/4]
    grid.lonq=grid.x_T_bounds[grid.jm/2,:]
    grid.latq=grid.y_T_bounds[:,grid.im/4]
    D=netCDF4.Dataset(grid_path+'topog_v4.nc').variables['depth'][:,:]
    grid.D=D
    grid.wet[grid.D<10.]=0.

    #lon,lat = np.meshgrid(grid.lonh,grid.lath)
    #i1,j1=near2d(lon, lat, lon1, lat1)
    #i2,j2=near2d(lon, lat, lon2, lat2)
    #i=[]; j=[]
    #i.append(i1);i.append(i2); i.sort()
    #j.append(j1);j.append(j2); j.sort()
    #region=grid.indexed_region(i=(i[0],i[1]),j=(j[0],j[1]))
    region=grid.geo_region(x=(lon1,lon2),y=(lat1,lat2)) 
    new_grid=grid.extract(region)

    #new_grid.lonh = new_grid.x_T[new_grid.jm/2,:]
    #new_grid.lonq = new_grid.x_T_bounds[new_grid.jm/2,:]
    #new_grid.lath = new_grid.y_T[:,new_grid.im/4]
    #new_grid.latq = new_grid.y_T_bounds[:,new_grid.im/4]

    
    ##section=grid.indexed_region(i=(ilat[0],ilat[-1]),j=(ilon[0],ilon[-1]))
    #section=grid.geo_region(x=(lonmin,lonmax),y=(latmin,latmax))
    #new_grid = grid.extract(section)        
    return new_grid,region

def get_depth(ilat,ilon):
    grid_path='/net2/mjh/ipynb/GIS_0125/'
    D=netCDF4.Dataset(grid_path+'topog_v4.nc').variables['depth'][ilat,ilon]
    return D

def near2d(x, y, x0, y0):
    """
    Find the indexes of the grid point that is
    nearest a chosen (x0, y0).
    Usage: line, col = near2d(x, y, x0, y0)
    """
    dx = np.abs(x - x0); dx = dx / dx.max()
    dy = np.abs(y - y0); dy = dy / dy.max()
    dn = dx + dy    
    fn = np.where(dn == dn.min())
    line = int(fn[0])
    col  = int(fn[1])
    return line, col

def get_area(var,dx,dy):
    [im,jm]=var.shape
    tot_area=0.;area=np.zeros((im,jm))
    area = np.ma.masked_where(np.ma.getmask(var), area) # applies the mask of m on x
    for i in range(im):
        for j in range(jm):
            if var[i,j] < 1.e20:
                area[i,j]=dx[i,j]*dy[i,j]
                tot_area=tot_area+(dx[i,j]*dy[i,j])
    return tot_area,area

def get_melt_rate(var,area,total_area):
    [jm,im]=var.shape
    melt_rate=np.zeros((jm,im))
    for i in range(im):
        for j in range(jm):
            melt_rate[j,i]=var[j,i]*area[j,i]/total_area

    return np.nansum(melt_rate)

def interm_pt(pnear, pk, pai, pbi, paj, pbj):
    ### FIND THE BEST INTERMEDIATE POINT ON A PATHWAY
    #           -----------------------------
    #   pnear   : vector of the position of the nearest point
    #   pk      : current working index
    #   pai, pbi: slope and original ordinate of x(y)
    #   paj, pbj: slope and original ordinate of y(x)
    #   pneari  : vector holding the position of intermediate point
    #           -----------------------------

    # 1 - Compute intermediate point

    # Determine whether we use y(x) or x(y):
    if (abs(paj) <= 1):
        # y(x)
        # possible intermediate point
        ylptmp1 = pnear[pk-1] + 1
        ylptmp2 = pnear[pk-1] + (paj/abs(paj))*1j
        # M is the candidate point:
        zxm = np.real(ylptmp1)
        zym = np.imag(ylptmp1)
        za0 = paj
        zb0 = pbj
        #
        za1 = -1./za0
        zb1 = zym-za1*zxm
        # P is the projection of M in the strait line
        zxp = -(zb1-zb0)/(za1-za0)
        zyp = za0*zxp+zb0
        # zd1 is the distance MP
        zd1 = (zxm-zxp) * (zxm-zxp) + (zym-zyp) * (zym-zyp)
        #
        # M is the candidate point:
        zxm = np.real(ylptmp2)
        zym = np.imag(ylptmp2)
        za1 = -1./za0
        zb1 = zym-za1*zxm
        # P is the projection of M in the strait line
        zxp = -(zb1-zb0)/(za1-za0)
        zyp = za0*zxp+zb0
        # zd1 is the distance MP
        zd2 = (zxm-zxp) * (zxm-zxp) + (zym-zyp) * (zym-zyp)
        #
        # choose the smallest (zd1,zd2)
        if (zd2 <= zd1):
            pneari = ylptmp2
        else:
            pneari = ylptmp1
        #       
    else:
        # x(y)
        ylptmp1 = pnear[pk-1] + (pai/abs(pai))
        ylptmp2 = pnear[pk-1] + 1*1j
        # M is the candidate point:
        zxm = np.real(ylptmp1)
        zym = np.imag(ylptmp1)
        za0 = pai
        zb0 = pbi
        #
        za1 = -1./za0
        zb1 = zxm-za1*zym
        # P is the projection of M in the strait line
        zyp = -(zb1-zb0)/(za1-za0)
        zxp = za0*zyp+zb0
        # zd1 is the distance MP
        zd1 = (zxm-zxp) * (zxm-zxp) + (zym-zyp) * (zym-zyp)
        #
        # M is the candidate point:
        zxm = np.real(ylptmp2)
        zym = np.imag(ylptmp2)
        za1 = -1./za0
        zb1 = zxm-za1*zym
        # P is the projection of M in the strait line
        zyp = -(zb1-zb0)/(za1-za0)
        zxp = za0*zyp+zb0
        # zd2 is the distance MP
        zd2 = (zxm-zxp) * (zxm-zxp) + (zym-zyp) * (zym-zyp)
        #
        # choose the smallest (zd1,zd2)
        if (zd2 <= zd1):
            pneari = ylptmp2
        else:
            pneari = ylptmp1
        
    return pneari

def section_transport(u, v, istart, iend, jstart, jend):
    """
    transpu, transpv = section_transport(u, v, grd,z_w, istart, iend, jstart, jend)
    compute the transport through the section defined between
    the point P1 (istart,jstart) and P2 (iend, jend).
    P1 and P2 are Arakawa-C psi points. z_w is the interface hight relative to mean sea level
    The transpot is positive right handside of the section.
    """


    # Find the nearest point between P1 (imin,jmin) and P2 (imax, jmax)
    # -----------------------------------------------------------------
    # Initialization
    i0=istart; j0=jstart; i1=iend;  j1=jend
    istart = float(istart); iend = float(iend)
    jstart = float(jstart); jend = float(jend)

    # Compute equation:  j = aj i + bj
    if istart != iend:
        aj = (jend - jstart ) / (iend - istart)
        bj = jstart - aj * istart
    else:
        aj=10000.
        bj=0.

    # Compute equation:  i = ai j + bi
    if jstart != jend:
        ai = (iend - istart ) / ( jend - jstart )
        bi = istart - ai * jstart
    else:
        ai=10000.
        bi=0.

    # Compute the integer pathway:
    # Chose the strait line with the smallest slope
    if (abs(aj) <=  1 ):
        # Here, the best line is y(x)
        print 'Here, the best line is y(x)'
        # If i1 < i0 swap points and remember it has been swapped
        if i1 <  i0:
            i  = i0 ; j  = j0
            i0 = i1 ; j0 = j1
            i1 = i  ; j1 = j
            norm = -1
        else:
            norm = 1

        if j1 >= j0:
            ist = 1; jst = 1
            norm_u = 1; norm_v = -1
        else:
            ist = 1; jst = 0
            norm_u = -1; norm_v = -1

        near = []
        # compute the nearest j point on the line crossing at i
        for i in range(i0,i1+1):
            j = aj*i + bj
            near.append(i + round(j)*1j)

    else:
        # Here, the best line is x(y)
        print 'Here, the best line is x(y)'
        # If j1 < j0 swap points and remember it has been swapped
        if j1 <  j0:
            i  = i0 ; j  = j0
            i0 = i1 ; j0 = j1
            i1 = i  ; j1 = j
            norm = -1
        else:
            norm = 1

        if i1 >= i0:
            ist = 1; jst = 1
            norm_u = 1; norm_v = -1
        else:
            ist = 0; jst = 1
            norm_u = 1; norm_v = 1

        near = []
        # compute the nearest i point on the line crossing at j
        for j in range(j0,j1+1):
            i = ai*j + bi
            near.append(round(i) + j*1j)

    # Look for intermediate points to be added
    # -------------------------------------------------------------

    inear = np.copy(near)

    n = len(near)
    nn=1

    for k in range(1,n):
        # distance between 2 neighbour points
        d = abs(inear[k] - inear[k-1])

        if ( d > 1 ):
            # intermediate points required if d>1
            neari = interm_pt(inear, k, ai, bi, aj, bj)
            near.insert(nn,neari)
            nn=nn+1

        nn=nn+1


    # Now extract the transport through a section
    # -------------------------------------------
    #set u and v to zero where u and v are masked for the sum
    for k in range(u.shape[0]):
        u[k,:] = np.where(u[k,:].mask, 0, u[k,:])
        v[k,:] = np.where(v[k,:].mask, 0, v[k,:])

    n = len(near)
    transpu = 0
    transpv = 0

    for l in range(0,n-1):
        ii = int(np.real(near[l])); jj = int(np.imag(near[l]))
        for k in range(0, u.shape[0]):
            if np.real(near[l]) == np.real(near[l+1]):
                trans = u[k, jj+jst, ii] * norm_u * norm
                transpu = transpu + trans

            elif np.imag(near[l]) == np.imag(near[l+1]):
                trans = v[k, jj, ii+ist] * norm_v * norm
                transpv = transpv + trans

    #plt.figure()
    #plt.plot(np.real(near),np.imag(near),'o')
    #plt.grid()
    #plt.xlim(np.min(np.real(near)-1),np.max(np.real(near)+1))
    #plt.ylim(np.min(np.imag(near)-1),np.max(np.imag(near)+1))
    #plt.show()

    return transpu, transpv
