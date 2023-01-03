import numpy as np
from astropy.io import fits
from mpdaf.obj import Spectrum,Cube,WCS,WaveCoord,iter_spe
import astropy.units as u
from numpy.polynomial import Polynomial as poly
import matplotlib.pyplot as plt


import os
import argparse 

parser = argparse.ArgumentParser()

parser=argparse.ArgumentParser()

parser.add_argument("-g", "--Galaxy", help = "Galaxy name")
parser.add_argument("-p", "--Path", help = "Path")
parser.add_argument("-z", "--Redshift", help = "Redshift")



args = parser.parse_args()

galaxy = args.Galaxy
path = args.Path
z = args.Redshift


os.chdir(path)
filename_emi = "{}_subtracted.fits".format(galaxy)
emi = Cube(filename_emi)

# we choose a random pixel just to get the wavelength-number relation

sp1 = emi[:,200,200]

# define redshift
oiii = 5008*z + 5008

# consider velocity dispersion for z1,z2 (+/-500 km/s)

z1,z2= [z-0.0016,z+0.0016]


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


