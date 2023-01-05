from astropy.io import fits
import numpy as np
import argparse 

parser = argparse.ArgumentParser()

parser=argparse.ArgumentParser()

parser.add_argument("-g", "--Galaxy", help = "Galaxy name")
parser.add_argument("-p", "--Path", help = "Path")


args = parser.parse_args()

galaxy = args.Galaxy
path = args.Path

path = path


for i in range(13):
    str_tmp = galaxy+'-{}_sli.fits'.format(i)
    hdu = fits.open(path+str_tmp)
    dt = hdu[1].data
    dt[dt==0] = np.nan
    hdu[1].data = dt
    str2_tmp = galaxy+'-{}_slinan.fits'.format(i)
    hdu.writeto(path+str2_tmp)
