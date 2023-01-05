from astropy.io import fits 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
import matplotlib as mpl
from photutils.aperture import CircularAperture,aperture_photometry
from scipy.stats import norm
from scipy.optimize import curve_fit
from numpy import asarray as ar,exp
import scipy.integrate as integrate
import os

import argparse 

parser = argparse.ArgumentParser()

parser=argparse.ArgumentParser()

parser.add_argument("-g", "--Galaxy", help = "Galaxy name")
parser.add_argument("-p", "--Path", help = "Path")
parser.add_argument("-z", "--Redshift", help = "Redshift", type=float)



args = parser.parse_args()

galaxy = args.Galaxy
path = args.Path
z = args.Redshift

os.chdir(path)


z1,z2= [z-0.0016,z+0.0016]


oiii_1 =  5008*z1 + 5008
oiii_2 = 5008*z2+5008

# we read the output of sextractor

spec = {}

for i in range(13):
    sli = 0+i
    cat_tmp = '{}-{}.cat'.format(galaxy,i)
    str_tmp = pd.read_csv(cat_tmp,delim_whitespace=True,comment='#',names=["x","y","ra","dec","mag","x2","y2","class"])
    nm_tmp = 'sli_{}'.format(i)
    spec[nm_tmp] = str_tmp
    

# matching since the criteria says we need a source to appear in at least three consecutive slices     
cat = {}

for m in range(0,11):
    
    cat_cen='sli_{}'.format(m+1)
    x,y = spec[cat_cen]['x'],spec[cat_cen]['y']
    mag_cen = spec[cat_cen]['mag']
    ra,dec = spec[cat_cen]['ra'],spec[cat_cen]['dec']
    
    cat_left='sli_{}'.format(m)
    x0,y0 = spec[cat_left]['x'],spec[cat_left]['y']
    mag_left = spec[cat_left]['mag']
    ra0,dec0 = spec[cat_left]['ra'],spec[cat_left]['dec']
    
    cat_right='sli_{}'.format(m+2)
    x1,y1 = spec[cat_right]['x'],spec[cat_right]['y']
    mag_right = spec[cat_right]['mag']
    ra1,dec1 = spec[cat_right]['ra'],spec[cat_right]['dec']
    
    dist = {}
    cat_t = []
    cat_nm = 'cat_{}'.format(m+1)
    
    for i in range(len(x)):
            dist0 = np.sqrt(np.power((x[i]-x0),2)+np.power((y[i]-y0),2))
            dist1 = np.sqrt(np.power((x[i]-x1),2)+np.power((y[i]-y1),2))
            if((np.any(dist0<0.5))and(np.any(dist1<0.5))):
                a = np.where(dist0<0.5)
                b = np.where(dist1<0.5)

                if((len(a)>1)or(len(b)>1)):
                    print('Error in',x[i],y[i])
                    
                else:
# we use the positions corresponding to the largest magnitude and assume thats the centroid

                    ar_cen = [x[i],y[i],ra[i],dec[i],mag_cen[i]]
                    ar_left = [x0[a[0][0]],y0[a[0][0]],ra0[a[0][0]],dec0[a[0][0]],mag_left[a[0][0]]]
                    ar_right = [x1[b[0][0]],y1[b[0][0]],ra1[b[0][0]],dec1[b[0][0]],mag_right[b[0][0]]]
            
                    mag = [ar_cen,ar_left,ar_right]
                    mags = np.array(mag)
                    which = np.where(max(mags[:,4]))
                    values = mags[which][0]
                
            
                    cat_t.append(values)
    cat[cat_nm] = cat_t
    
# we check for repetitions
cati = np.array(cat['cat_1'])
x0 = cati[:,0]
y0 = cati[:,1]

cat_tmp = cat['cat_1']

for m in range(2,11):
    cat_nm = 'cat_{}'.format(m)
    for i in range(len(cat[cat_nm])):
        x = cat[cat_nm][i][0]
        y = cat[cat_nm][i][1]
        
        dist = 0.5
        if(np.all(np.sqrt(np.power((x-x0),2)+np.power((y-y0),2))>dist)):
            values = cat[cat_nm][i]
            cat_tmp.append(values)
            
cat_fin = np.array(cat_tmp)


# we perform aperture photometry
x,y = cat_fin[:,0],cat_fin[:,1]
ra,dec = cat_fin[:,2],cat_fin[:,3]

f = len(x)

positions = []
pos_fin = []
for i in range(len(x)):
    positions.append([x[i],y[i]])
    pos_fin.append([ra[i],dec[i],x[i],y[i]])

table = {}
for i in range(13):
    sliu = '{}-{}_sliun.fits'.format(galaxy,i)
    tmp = 'phot_{}'.format(i)
    aperture = CircularAperture(positions,r=3.)
    hdu = fits.open(sliu)
    data = hdu[1].data
    table[tmp] = aperture_photometry(data,aperture)

# we fit gaussians to make sure there is a signal
aux = []
tot = len(table[tmp]['xcenter'].value)

for k in range(2,3):
    for i in range(13):
        tmp = 'phot_{}'.format(i)
        aux.append(table[tmp][k]['aperture_sum'])
wave = np.arange(oiii_1,oiii_2-1.25,1.25).tolist()


sp = np.array(aux[0:12])
wv = np.array(wave[0:12])

mean,std = norm.fit(wv)
x = wv
y = sp


def gaus(x,a,me,st,g):
    return a*exp(-(x-me)**2/(2*st**2))+g

flag= np.ones(f)
det = np.zeros(f)
flux = np.zeros(f)
err = np.zeros(f)
cen = np.zeros(f)
tmp = []
spec= []
vrad = np.zeros(f)
oiii = (z+1)*5008

g = 0
for k in range(len(flag)):
    spec = []
    for i in range(13):
        tmp = 'phot_{}'.format(i)
        spec.append(table[tmp][k]['aperture_sum'])

    wave = np.arange(oiii_1,oiii_2-1.25,1.25).tolist()
    
    sp = np.array(spec[2:10])
    wv = np.array(wave[2:10])
    mean,std = norm.fit(wv)
    x = wv
    y = sp

    try:
        popt,pcov = curve_fit(gaus,x,y,p0=[max(sp),mean,std,g],bounds=((0,5000,-np.inf,-np.inf),(10000,5500,np.inf,np.inf)))
        det[k] = np.linalg.det(pcov)
        
        if((det[k]==0)or(np.log10(det[k])>5)or(popt[3]>0)):
            flag[k] = 1
        else:
            tmp = integrate.quad(gaus,oiii_1,oiii_2,args=(popt[0],popt[1],popt[2],0))
            flux[k] = tmp[0]
            err[k] = tmp[1]
            flag[k] = 0
            cen[k] = popt[1]
            vrad[k] = 300000*(popt[1]-oiii)/oiii
    except:        
        flag[k] = 1
        det[k] = np.nan

        
final = []
for i in range(len(flag)):
    if(flag[i]==0):
        final.append([pos_fin[i][0],pos_fin[i][1],pos_fin[i][2],pos_fin[i][3],flux[i],cen[i]])
    
fin = np.array(final)


fin = pd.DataFrame(fin,columns=['ra','dec','x','y','flux','cen'])

ar = []
for i in range(len(fin['x'])):
    x = fin['x'][i]
    y = fin['y'][i]
    x0 = fin['x']
    y0 = fin['y']
    if(np.any((np.sqrt(np.power((x-x0),2)+np.power((y-y0),2)))<0.5)):
        ind = np.where((np.sqrt(np.power((x-x0),2)+np.power((y-y0),2)))<0.5)
        for j in range(len(ind[0])):
            if(ind[0][j]!=i):
                if(np.all(ind[0][j]!=ar)):
                    ar.append(i)   
            
fin = fin.drop(ar)
fin = fin.reset_index(drop=True)
fin = fin.sort_values(by=['x','y'])
fin = fin.reset_index(drop=True)
 

catalogue = "./{}_pne.cat".format(galaxy)
np.savetxt(catalogue,fin)


ra = cat_fin[:,2]
dec = cat_fin[:,3]
cat_vrad = np.transpose([ra,dec,vrad])
catv = "./{}_vrad_wcs.cat".format(galaxy)
np.savetxt(catv,cat_vrad)