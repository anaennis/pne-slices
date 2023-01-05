import os
from astropy.io import fits
import argparse 

parser = argparse.ArgumentParser()

parser=argparse.ArgumentParser()

parser.add_argument("-g", "--Galaxy", help = "Galaxy name")
parser.add_argument("-p", "--Path", help = "Path")



args = parser.parse_args()

galaxy = args.Galaxy
path = args.Path

os.chdir(path)

cube = "./../{}.fits".format(galaxy) 
image = fits.open(cube) #read MUSE cube or DATACUBE_FINAL
#extract useful values from cube header
CRPIX1 = image[1].header['CRPIX3']
CRDEL1 = image[1].header['CD3_3']
CRVAL1 = image[1].header['CRVAL3']


emi = "./../{}emi.fits".format(galaxy) 
image2 = fits.open(cube) #read MUSE cube or DATACUBE_FINAL
#extract useful values from cube header
CRPIX2 = image2[1].header['CRPIX3']
CRDEL2 = image2[1].header['CD3_3']
CRVAL2 = image2[1].header['CRVAL3']


important_data = [cube,CRPIX1,CRDEL1,CRVAL1,emi,CRPIX2,CRDEL2,CRVAL2]
dat = "./{}_dat.cat".format(galaxy)

with open(dat, 'w') as f:
    for line in important_data:
        f.write(f"{line}\n")