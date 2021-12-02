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

styling = {
        'text.usetex' : True,
        'lines.linewidth' : 2,
        'font.family' : 'serif',
        'axes.prop_cycle' : cycler('color', ['#EE6666', '#3388BB', '#9988DD',
                 '#EECC55', '#88BB44', '#FFBBBB']),
        'axes.labelsize' : 20,
        'legend.fontsize' : 15,
        "xtick.labelsize": 15,
        "ytick.labelsize": 15,
        'ytick.major.size' : 5,
        'ytick.minor.size' : 4,
        'xtick.major.size' : 5,
        'xtick.minor.size' : 4,
        'xtick.direction' : 'inout',
        'ytick.direction' : 'inout',  
}

mpl.rcParams.update(styling)

name = 'ngc1387'
z = 0.0043
oiii_1 = 5021.5216 
oiii_2 = 5038.048

# we read the output of sextractor

spec = {}

for i in range(15):
    sli = 0+i
    cat_tmp = 'slice-{}.cat'.format(i)
    str_tmp = pd.read_csv(cat_tmp,delim_whitespace=True,comment='#',names=["x","y","ra","dec","mag","x2","y2","class"])
    nm_tmp = 'sli_{}'.format(i)
    spec[nm_tmp] = str_tmp
    

# matching since the criteria says we need a source to appear in at least three consecutive slices     
cat = {}

for m in range(0,13):
    
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

for m in range(2,13):
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

positions = []
pos_fin = []
for i in range(len(x)):
    positions.append([x[i],y[i]])
    pos_fin.append([ra[i],dec[i],x[i],y[i]])

table = {}
for i in range(13):
    sliu = '{}-{}_sliun.fits'.format(name,i)
    tmp = 'phot_{}'.format(i)
    aperture = CircularAperture(positions,r=3.)
    hdu = fits.open(sliu)
    data = hdu[1].data
    table[tmp] = aperture_photometry(data,aperture)

# we fit gaussians to make sure there is a signal
aux = []
tot = len(table[tmp]['xcenter'].value)

for k in range(25,26):
    for i in range(13):
        tmp = 'phot_{}'.format(i)
        aux.append(table[tmp][k]['aperture_sum'])
wave = np.arange(oiii_1,oiii_2-1.25,1.25).tolist()


sp = np.array(aux[2:12])
wv = np.array(wave[2:12])

mean,std = norm.fit(wv)
x = wv
y = sp


def gaus(x,a,me,st,g):
    return a*exp(-(x-me)**2/(2*st**2))+g

flag= np.ones(67)
det = np.zeros(67)
flux = np.zeros(67)
err = np.zeros(67)
cen = np.zeros(67)
tmp = []
spec= []
vrad = np.zeros(67)
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
        
        if((det[k]==0)or(np.log10(det[k])>3.5)or(popt[3]>0)):
            flag[k] = 1
        else:
            tmp = integrate.quad(gaus,5020,5040,args=(popt[0],popt[1],popt[2],0))
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

catalogue = "./{}_pne.cat".format(name)
np.savetxt(catalogue,fin)


ra = cat_fin[:,2]
dec = cat_fin[:,3]
cat_vrad = np.transpose([ra,dec,vrad])
catv = "./{}_vrad_wcs.cat".format(name)
np.savetxt(catv,cat_vrad)