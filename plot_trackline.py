from midas.rectgrid import *
from midas.wright_eos import *
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import sys
sys.path.append('/net/gmm/Libs/Python/mycache/')
from cache import *
sys.path.append('/net/gmm/Libs/Python/trackline/')
from trackline import *
import argparse

global nclick

nclick=0

def MFlist(d):

    import glob
    
    files=[]
    pre='gfdl.ncrc2-default-prod/history'
    listing= d['root']+d['exp']+'/'+d['name']+'/'+pre+'/'+d['suffix']
    listing = glob.glob(listing)

    for filename in listing:
        files.append(filename)
    
    files=sorted(files)
    return files

#### Begin User Input

parser = argparse.ArgumentParser()
parser.add_argument('-n',type=int,help='time slice',default=0)
parser.add_argument('--exp',type=str,help='experiment',default=None)
parser.add_argument('--trackname',type=str,help='previously saved track',default=None)
parser.add_argument('--region',type=str,help='FR|ROSS|AMERY|GIBRALTAR|LAB|NWATL|SEPAC|GIS',default=None)
parser.add_argument('--field',type=str,help='temp|salt|age|Kd_work',default=None)
parser.add_argument('--vmin',type=float,help='min val for plotting',default=None)
parser.add_argument('--vmax',type=float,help='max val for plotting',default=None)
parser.add_argument('--ignore',type=float,help='missing value',default=None)
parser.add_argument('--ignore_lt',type=float,help='missing value',default=None)
parser.add_argument('--ignore_gt',type=float,help='missing value',default=None)
parser.add_argument('--zlim',type=float,help='vertical section depth (m)',default=None)
parser.add_argument('--z0',type=float,help='vertical section depth (m)',default=0.0)
parser.add_argument('--savefig',type=bool,help='save image',default=False)

args=parser.parse_args()

def create_trackline(P,S,field):
    dict={}
    dict['iind']=[]
    dict['jind']=[]
    dict['dates']=[]
    dict['lon']=[]
    dict['lat']=[]
    dict[field]=[]
    dict['z']=[]
    dict['zi']=[]

    for p in P:
        dates=S.var_dict[field]['dates'].copy()
        dates=instance_to_datetime(dates)
        dict['jind'].append(p[0])
        dict['iind'].append(p[1])
        dict['dates'].append(dates)
        dict['lon'].append(grid.x_T[p[0],p[1]])
        dict['lat'].append(grid.y_T[p[0],p[1]])
        dict[field].append(sq(vars(S)[field][:,:,p[0],p[1]]))
        dict['z'].append(sq(S.var_dict[field]['z'][:,:,p[0],p[1]]))
        dict['zi'].append(sq(S.var_dict[field]['z_interfaces'][:,:,p[0],p[1]]))
        
    return dict

def update_trackline(d,P,S,field):
    d[field]=[]

    for p in P:
        d[field].append(sq(vars(S)[field][:,:,p[0],p[1]]))
 

    return field

def read_trackline(fnam):

    import csv

    P=[]
    f=open(fnam,'r')
    a=csv.reader(f)
    for line in a:
        P.append(map(int,line))

    return P
# Load Grid and topography

sgrid=supergrid(file='/net3/mjh/models/GIS_0125/0125gridGeneration/ocean_hgrid.nc',cyclic_x=True,tripolar_n=True)
grid=quadmesh(supergrid=sgrid)
grid.D=netCDF4.Dataset('/net3/mjh/projects/GIS_0125/topog_v4.nc').variables['depth'][:]
grid.wet[grid.D==0.]=0.
grid.lath=grid.y_T[:,grid.im/4]
grid.latq=grid.y_T_bounds[:,grid.im/4]
IS = nc.Dataset('/net3/mjh/projects/GIS_0125/ice_shelf_v4.nc').\
variables['thick'][:]
IS[IS>0.]=1.


# Define regions

FR=grid.indexed_region(i=(1720,2410),j=(0,515))
ROSS=grid.indexed_region(i=(760,1320),j=(0,360))
AMERY=grid.indexed_region(i=(0,260),j=(225,445))
GIBRALTAR=grid.indexed_region(i=(2300,2380),j=(1365,1440))
NWATL=grid.geo_region(x=(-85,-10),y=(30,70))
BROKE=grid.geo_region(x=(40,75),y=(-74,-50))
SEPAC=grid.geo_region(x=(-120,-67),y=(-30,5))
BOB=grid.geo_region(x=(-285,-260),y=(5,25))
GIS=grid.indexed_region(i=(2040,2480),j=(1680,1960))


region_names=['FR','ROSS','AMERY','GIBRALTAR','NWATL','SEPAC','BOB','GIS']
regions=[FR,ROSS,AMERY,GIBRALTAR,NWATL,SEPAC,BOB,GIS]

# Get Ice Shelf thickness for regions (np.array)

FR_IS = nc.Dataset('/net3/mjh/models/GIS_0125/0125gridGeneration/ice_shelf/ice_shelf_v1.nc').\
variables['thick'][FR['y'],FR['x_read']]
ROSS_IS = nc.Dataset('/net3/mjh/models/GIS_0125/0125gridGeneration/ice_shelf/ice_shelf_v1.nc').\
variables['thick'][ROSS['y'],ROSS['x_read']]
AMERY_IS = nc.Dataset('/net3/mjh/models/GIS_0125/0125gridGeneration/ice_shelf/ice_shelf_v1.nc').\
variables['thick'][AMERY['y'],AMERY['x_read']]

reg_is = None

if args.region == 'FR':
    GEO_REGION=FR
    reg=grid.extract(FR)
    reg_is=FR_IS
elif args.region == 'ROSS':
    GEO_REGION=ROSS    
    reg=grid.extract(ROSS)
    reg_is=ROSS_IS
elif args.region == 'AMERY':
    GEO_REGION=AMERY
    reg=grid.extract(AMERY)
    reg_is=AMERY_IS
elif args.region == 'GIBRALTAR':
    GEO_REGION=GIBRALTAR
    reg=grid.extract(GIBRALTAR)
elif args.region == 'NWATL':
    GEO_REGION=NWATL
    reg=grid.extract(NWATL)    
elif args.region == 'SEPAC':
    GEO_REGION=SEPAC
    reg=grid.extract(SEPAC)    
elif args.region == 'BOB':
    GEO_REGION=BOB
    reg=grid.extract(BOB)    
elif args.region == 'GIS':
    GEO_REGION=GIS
    reg=grid.extract(GIS)    
elif args.region == 'BROKE':
    GEO_REGION=BROKE
    reg=grid.extract(BROKE)    

if reg_is is not None:
    reg_is[reg_is>0]=1

GEO_REGION=ROSS
reg=grid.extract(ROSS)
reg_is=ROSS_IS
mask=np.zeros(reg.wet.shape)
mask=mask.astype(int)
    
if args.trackname is None:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.pcolormesh(np.ma.masked_where(mask==1,mask+1),alpha=0.25)
    ax.contour(reg.D,[0.,100,200,500,1000,2000],colors='k')
    if reg_is is not None:
        ax.contour(reg_is,[0.5,0.5],colors='purple')
    xdata=[]
    ydata=[]

    def onclick(event):

        global nclick
        
        print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
            event.button, event.x, event.y, event.xdata, event.ydata)
        ident=str(nclick)
        nclick=nclick+1
        ax.text(event.xdata, event.ydata,ident, fontsize=8)
        xdata.append(int(event.xdata))
        ydata.append(int(event.ydata))

        fig.canvas.draw()
        
    cid = fig.canvas.mpl_connect('button_press_event', onclick)

    plt.show()

    fnam='.track'
    f=open(fnam,'w')

    for n in np.arange(0,len(xdata)):
        a=str(xdata[n])+','+str(ydata[n])+'\n'
        f.write(a)
    
    f.close()
else:
    fnam='.track_'+args.trackname
    P=read_trackline(fnam)
    xdata=[];ydata=[]
    for p in P:
        xdata.append(p[0])
        ydata.append(p[1])

        
mp=np.zeros(mask.shape)
path=[]
for n in np.arange(0,len(xdata)-1):
    p,m=node_path(xdata[n],ydata[n],xdata[n+1],ydata[n+1],mask,dirs=4,metric='Manhattan')
    for node in p:
        path.append(node)
    mp=mp+np.array(m)
    mp[mp==1]=0


fig = plt.figure()
ax = fig.add_subplot(111)

ax.pcolormesh(reg.x_T,reg.y_T,np.ma.masked_where(mp==0,mp),cmap=plt.cm.flag,alpha=0.5)
if reg_is is not None:
    ax.contour(reg.x_T,reg.y_T,reg_is,[0.5,0.5],colors='b')
ax.contour(reg.x_T,reg.y_T,reg.D,[0.,100,200,500,1000,2000],colors='k')

d={}
d['root']='/ptmp/Matthew.Harrison/archive/Matthew.Harrison/fre/ulm_prerelease/'
#d['root']='/ptmp/Matthew.Harrison/archive/Matthew.Harrison/fre/tikal_201407/'
d['exp']='GIS_0125_LM3_SIS'
#d['name']=args.exp
d['name']='GIS_0125_LM3_SIS_shelfthermo4/'
d['suffix']='*/*.ocean_month_snap.nc'
fnames=MFlist(d)


T=state(MFpath=fnames,fields=['temp'],grid=grid,geo_region=GEO_REGION,interfaces='e',time_indices=[args.n],verbose=False)

args.field='temp'

if args.field in ['age','Kd_work']:
    d['suffix']='*/*.ocean_month.nc'
    fnames=MFlist(d)
    F=state(MFpath=fnames,fields=[args.field],grid=grid,geo_region=GEO_REGION,time_indices=[args.n],verbose=False)
elif args.field in ['temp','salt','u','v']:
    d['suffix']='*/*.ocean_month_snap.nc'
    fnames=MFlist(d)
    print fnames
    F=state(MFpath=fnames,fields=[args.field],grid=grid,geo_region=GEO_REGION,time_indices=[args.n],verbose=False)
else:
    print args.field, ' not available'
    raise()

print 'args.n',args.n

d=create_trackline(path,T,'temp')
update_trackline(d,path,F,args.field)

zout=np.array(d[args.field]).T

if args.ignore is not None:
    zout=np.ma.masked_where(zout==args.ignore,zout)

if args.ignore_lt is not None:
    zout=np.ma.masked_where(zout<args.ignore_lt,zout)

if args.ignore_gt is not None:
    zout=np.ma.masked_where(zout>args.ignore_gt,zout)    

x=np.arange(0,len(d['lon']))
#x=0.5*(np.array(d['lon'][0:-1])+np.array(d['lon'][1:]))
#x=np.concatenate(([d['lon'][0]-0.125],x))
Z=np.array(d['zi'][:]).T

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)

X,z=np.meshgrid(x,np.arange(0,64))

if args.vmin is not None:
    cf=ax2.pcolormesh(X,Z,zout,vmin=args.vmin,vmax=args.vmax,cmap=plt.cm.jet)
else:
    cf=ax2.pcolormesh(X,Z,zout,cmap=plt.cm.jet)    

if args.zlim is not None:
    plt.ylim(-args.zlim,args.z0)


zbot=Z[-1,:]
hml=Z[2,:]
hbl=Z[4,:]
hml=np.ma.masked_where(zbot>-10.,hml)
hbl=np.ma.masked_where(zbot>-10.,hbl)
zbot=np.ma.masked_where(zbot>-10.,zbot)

plt.plot(x,zbot,color='k',linewidth=2.0)
plt.plot(x,Z[0,:],color='k',linewidth=1.0)
plt.plot(x,hml,color='b',linewidth=2.0)
plt.plot(x,hbl,color='g',linewidth=2.0)

plt.xlim(0,X.shape[1])
plt.grid()
plt.colorbar(cf)
#tit=args.exp+' '+args.field+' '+datetime.datetime.ctime(d['dates'][0][0])
#plt.title(tit,fontsize=10)
if args.savefig:
    track=args.trackname
    if track is None: track='track'
    fnam=track+'_'+args.field+'_'+'%04d'%args.n +'.png'
    plt.savefig(fnam)
else:
    plt.show()

