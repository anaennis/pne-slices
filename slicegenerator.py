import numpy as np
from astropy.io import fits
from mpdaf.obj import Spectrum,Cube,WCS,WaveCoord,iter_spe
import astropy.units as u
from numpy.polynomial import Polynomial as poly
import matplotlib.pyplot as plt

from cycler import cycler
import matplotlib as mpl

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

name = "ngc1387"

filename = "./{}.fits".format(name) 
filename_cont = "./{}cont.fits".format(name) 
filename_im = "{}_ima.fits".format(name)

cube = Cube(filename)
cont = Cube(filename_cont)
emi = cube - cont

# we choose a random pixel just to get the wavelength-number relation

sp1 = cube[:,200,200]

# define redshift
z = 0.0043
oiii = 5008*z + 5008

# consider velocity dispersion for z1,z2
z1,z2= [0.0027,0.006]


oiii_1 =  5008*z1 + 5008
oiii_2 = 5008*z2+5008
ka,kb = sp1.wave.pixel([oiii_1,oiii_2],nearest=True)


spec = {}
spec2 = {}
for i in range(13):
    sli = ka+i
    tmp = emi[(sli-1):(sli+2),:,:]
    im_tmp = tmp.sum(axis=0)
    str_tmp = 'slice_{}'.format(i)
    spec[str_tmp] = im_tmp
    file_out = '{}-{}_sli.fits'.format(name,i)
    spec[str_tmp].write(file_out)
    
# we create unbinned slices to obtain the flux for a single slice later on    
    im_tmp2 = emi[sli,:,:]
    str_tmp2 = 'slice_unbin_{}'.format(i)
    spec2[str_tmp2] = im_tmp2
    file_out2 = '{}-{}_sliun.fits'.format(name,i)
    spec2[str_tmp2].write(file_out2)
    
    
filtro_oiii = emi[ka:kb,:,:]
im_oiii = filtro_oiii.sum(axis=0)
filename_out = "{}_ox.fits".format(name) 

im_oiii.write(filename_out)


