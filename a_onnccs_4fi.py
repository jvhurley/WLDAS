import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.path import Path


### data files
data1=('/discover/nobackup/projects/wldas/pub/wldas_domain/mthly_recharge/mthly_gw_rcg_bigdm1979-1989.bin')
#f=open(data1,"rb")

data2=('/discover/nobackup/projects/wldas/pub/wldas_domain/mthly_recharge/mthly_gw_rcg_bigdm1990-1999.bin')
data3=('/discover/nobackup/projects/wldas/pub/wldas_domain/mthly_recharge/mthly_gw_rcg_bigdm2000-2009.bin')
data4=('/discover/nobackup/projects/wldas/pub/wldas_domain/mthly_recharge/mthly_gw_rcg_bigdm2010-2017.bin')

### extent of wldas data
ndy=2787
ndx=3591
strt=0
step=(ndy*ndx)
edge_s=25.065
edge_n=52.925
edge_w=-124.925
edge_e=-89.025

### full lon,lat
lat_g=np.linspace(edge_s,edge_n,ndy)
lon_g=np.linspace(edge_w,edge_e,ndx)

### user-defined area (currently CA and NV)
gwb_enca=43.
gwb_esca=30.
gwb_ewca=lon_g[0]
gwb_eeca=-113.
gwb_gencai=(np.abs(lat_g-gwb_enca)).argmin()
gwb_gescai=(np.abs(lat_g-gwb_esca)).argmin()
gwb_gewcai=(np.abs(lon_g-gwb_ewca)).argmin()
gwb_geecai=(np.abs(lon_g-gwb_eeca)).argmin()

### user area lon,lat
lon_gca=lon_g[gwb_gewcai:gwb_geecai]
lat_gca=lat_g[gwb_gescai:gwb_gencai]

### number of years in data files
yinfs=np.array(([11,10,10,8]))

### choose data file to use, used below
def datanamfunc(fileno):
    if fileno==0:
        dafi=data1
    if fileno==1:
        dafi=data2
    if fileno==2:
        dafi=data3
    if fileno==3:
        dafi=data4
    return(dafi)

### read in the data and build the 4d array of year,month,lon,lat
yyc=0
for i in np.arange(0,yinfs.shape[0]):
    yinf=yinfs[i]
    data=datanamfunc(i)
    f=open(data,"rb")
    for yy in np.arange(0,yinf):
        for mm in np.arange(0,12):
            tgd=[]
            if yy==0 and mm==0:
                f.seek(strt,os.SEEK_SET)
                tgd=np.fromfile(f,dtype='<f',count=step,sep="")
                gw1=tgd.reshape(ndy,ndx)
                gw2=gw1[gwb_gescai:gwb_gencai,gwb_gewcai:gwb_geecai]
                if i==0:
                    gwmosy=np.zeros((np.nansum(yinfs),12,gw2.shape[0],gw2.shape[1]))
                gwmosy[yyc,mm,:,:]=gw2
                #gwmost=np.zeros((12,gw1.shape[0],gw1.shape[1]))
            else:
                f.seek(strt,os.SEEK_CUR)
                tgd=np.fromfile(f,dtype='<f',count=step,sep="")
                gw1=tgd.reshape(ndy,ndx)
                gw2=gw1[gwb_gescai:gwb_gencai,gwb_gewcai:gwb_geecai]
                gwmosy[yyc,mm,:,:]=gw2
                #gwmost[mm,:,:]=gwmost[mm,:,:]+gw1
            print(mm,np.nanmean(gw2))
        yyc=yyc+1
        print(yy)
    print(data)
    f.close()

### switch the no-value values to nans and take a mean across the years-->climatology
gwmosy[gwmosy<-1000]=np.nan
gwmosy_yt=np.nanmean(gwmosy,axis=0)

#### for work on discover with full time-series
#np.savez('/discover/nobackup/projects/wldas/jvhurley/CAgwrc_7917.npz',gwmosy=gwmosy,lon_gca=lon_gca,lat_gca=lat_gca)
#### for work on LOCAL with climo
#np.savez('/discover/nobackup/projects/wldas/jvhurley/CAgwrc_7917_yt.npz',gwmosy_yt=gwmosy_yt,lon_gca=lon_gca,lat_gca=lat_gca)






























####### PROTOTYPE FOR BASINSERIES 
acres2sqmet_div=0.00024711
npzf=('/discover/nobackup/projects/wldas/jvhurley/gwb_borval_borspr.npz')
gwb_bv=np.load(npzf)['gwb_borval_borspr']
gwb_bv_path=Path(gwb_bv)
gwb_bv_area_acres=62749.2
gwb_bv_area_sqmet=gwb_bv_area_acres/acres2sqmet_div

def basinseries(lon_big,lat_big,gwb,gwmy,gwbp):  ## lon,lat,gwbasin,gwrecharge(yy,mm,x,y),gwbasin-path
    ##### extent of gwbasin
    gwb_en=np.max(gwb[:,1])
    gwb_es=np.min(gwb[:,1])
    gwb_ew=np.min(gwb[:,0])
    gwb_ee=np.max(gwb[:,0])
    gwb_geni=(np.abs(lat_big-gwb_en)).argmin()
    gwb_gesi=(np.abs(lat_big-gwb_es)).argmin()
    gwb_gewi=(np.abs(lon_big-gwb_ew)).argmin()
    gwb_geei=(np.abs(lon_big-gwb_ee)).argmin()
    lon_gwb=lon_big[gwb_gewi:gwb_geei]
    lat_gwb=lat_big[gwb_gesi:gwb_geni]
    basints=np.zeros((gwmy.shape[0],gwmy.shape[1]))
    for aa in np.arange(gwmy.shape[0]):
        for bb in np.arange(gwmy.shape[1]):
            gwmos_tmp=gwmy[aa,bb,gwb_gesi:gwb_geni,gwb_gewi:gwb_geei]
            gwr_ab_xyz=np.zeros((gwmos_tmp.shape[0]*gwmos_tmp.shape[1],3))
            gwcr=0
            for i in np.arange(0,gwmos_tmp.shape[0]):
                for j in np.arange(0,gwmos_tmp.shape[1]):
                    gwr_ab_xyz[gwcr,0]=lon_gwb[j]
                    gwr_ab_xyz[gwcr,1]=lat_gwb[i]
                    gwr_ab_xyz[gwcr,2]=gwmos_tmp[i,j]
                    gwcr=gwcr+1
            if (aa==0) and (bb==0):
                gwrii=np.array(np.transpose(np.where(gwbp.contains_points(gwr_ab_xyz[:,0:2])==1)))
            basints[aa,bb]=np.nanmean(gwr_ab_xyz[gwrii,2],axis=0)
        print(aa)
    return(basints)

gwbr=basinseries(lon_gca,lat_gca,gwb_bv,gwmosy,gwb_bv_path)

cubmet_2_af_mult=0.000810714                       ## multiplier to convert cubic meters to acre feet
gwbr_af=((gwbr/1000.)*gwb_bv_area_sqmet)*cubmet_2_af_mult









####################### LOCAL LOCAL LOCAL LOCAL #####################  LOCAL  ############
#from matplotlib.path import Path
npzfa=('C:/Users/JHurley/Mangoes/nasa/CAgwrc_7917.npz')
gwmosy=np.load(npzfa)['gwmosy']
lon=np.load(npzfa)['lon_gca']
lat=np.load(npzfa)['lat_gca']
gwmosy_yt=np.nanmean(gwmosy,axis=0)
gwmosy_jan=gwmosy_yt[0,:,:]
xx,yy=np.meshgrid(lon,lat)

npzf=('C:/Users/JHurley/Mangoes/nasa/gwb_borval_borspr.npz')
gwb=np.load(npzf)['gwb_borval_borspr']
gwb_path=Path(gwb)

















