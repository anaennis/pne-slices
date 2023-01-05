#astropy-photutils-related functions
from photutils import CircularAperture, CircularAnnulus
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.stats import sigma_clipped_stats
from astropy.visualization import simple_norm
from photutils import aperture_photometry
from photutils import DAOStarFinder
from astropy.table import Table
from astropy.io import fits

#usual python imports
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import warnings
import glob
import os

import argparse 

parser = argparse.ArgumentParser()

parser=argparse.ArgumentParser()

parser.add_argument("-g", "--Galaxy", help = "Galaxy name")
parser.add_argument("-p", "--Path", help = "Path")



args = parser.parse_args()

galaxy = args.Galaxy
path = args.Path

os.chdir(path)
#function to generate an arbitrary wl range of MUSE. Full coverage: (3721)*1.25+4699.625
def wl(N, crpix, crdelt, crval):
    """
    N: numero del slice
    crpix: reference slice de inicio
    crdelt: delta lambda por cada slice
    crval: lambda de referencia/inicio
    """
    return (N+1-crpix)*crdelt + crval


catalogue = "{}_pne.cat".format(galaxy)
sources = pd.read_csv(catalogue,delim_whitespace=True,names=['ra','dec','xcentroid','ycentroid','flux','cen'])
pos = sources[['xcentroid', 'ycentroid']].values #pos contains only (x,y) positions 
apertures = CircularAperture(pos, r=5) #create the apertures for the pos pixel locations within a radius r
annulus_apertures = CircularAnnulus(pos, r_in=5, r_out=8) #create an annulus to be used as background 
annulus_masks = annulus_apertures.to_mask(method='center')

cube = "./../{}.fits".format(galaxy) 
image = fits.open(cube) #read MUSE cube or DATACUBE_FINAL
#extract useful values from cube header
CRPIX1 = image[1].header['CRPIX3']
CRDEL1 = image[1].header['CD3_3']
CRVAL1 = image[1].header['CRVAL3']


# This routine was adapted from the Astropy-photutils package.
# Performs aperture extractions + background substraction (pre-defined annulus within r_in and r_out from 2 cells above)
# https://photutils.readthedocs.io/en/stable/aperture.html

spec = pd.DataFrame()
slices_to_extract = np.arange(image[1].data.shape[0]) #change the number of slices to be extracted here
for n_slice in slices_to_extract: #receives the n-th slice index to be extracted
    print('Slice %d/%d' %(n_slice, slices_to_extract.size-1), end='\r')
    
    #extract n-th slice data/flux 
    data = image[1].data[n_slice,:,:]
    data = np.nan_to_num(data)
    data[data < 0] = 0
    #extract n-th slice error data
    err = image[2].data[n_slice,:,:]
    err = np.nan_to_num(err)
    err[err < 0] = 0
    
    bkg_median = []
    bkg_error = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(data)
        annulus_data_1d = annulus_data[mask.data > 0]
        mean_sigclip, median_sigclip, std_sigclip = sigma_clipped_stats(annulus_data_1d)
        sqN = np.sqrt(np.sum(mask.data > 0))
        bkg_median.append(median_sigclip)
        bkg_error.append(1.253*std_sigclip/sqN)

    bkg_median = np.array(bkg_median)
    bkg_error = np.array(bkg_error)
    phot = aperture_photometry(data, apertures, error=err**0.5)
    phot['annulus_median'] = bkg_median
    phot['aper_bkg'] = bkg_median * apertures.area
    phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
    phot['total_error'] = np.sqrt((phot['aperture_sum_err'])**2.+(bkg_error)**2.)
    phot = pd.DataFrame(phot.as_array())
    
    if n_slice == 0:
        labels = ['id', 'xcenter', 'ycenter', 'aper_sum_bkgsub', 'total_error']
        spec = phot[labels]
    else:
        labels = ['id', 'aper_sum_bkgsub', 'total_error']
        spec = pd.merge(spec, phot[labels], on='id')
        
nm_fin = '{}_aper_phot.csv'.format(galaxy)
spec.to_csv(nm_fin, index=False)
#each coma separated column after the stars id, xcenter, and ycenter contains the flux and error of each slice