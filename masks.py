from astropy.io import fits
import numpy as np

# open any cube
maskvieja = fits.open('ngc1387_oiii-Copy1.fits')


# create array with same size as image
# arrays are defined as (y,x)
array = np.zeros((318,318))

# define masks

for y in range(295,311):
    for x in range(119,138):
        array[y][x] = 1

        
for y in range(320,329):
    for x in range(234,244):
        array[y][x] = 1


# set up format
prim = maskvieja[0].data
primary_hdu = fits.PrimaryHDU(prim)
image_hdu = fits.ImageHDU(array)
hdul = fits.HDUList([primary_hdu,image_hdu])

hdul.writeto('mask5128.fits')


