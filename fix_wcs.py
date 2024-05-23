from astropy.wcs import WCS
from astropy.io import fits
import argparse
import os
from mpdaf.obj import Cube

parser = argparse.ArgumentParser()

parser.add_argument("-g","--Galaxy",help="Galaxy Name")
parser.add_argument("-p","--Path",help="Path to the slices folder")

args= parser.parse_args()

galaxy = args.Galaxy
path = args.Path

os.chdir(path)
nm = 'stitched.fits'.format(galaxy)

hd = fits.open(nm)


len1 = hd[1].header['NAXIS2']
len2 = hd[1].header['NAXIS1']

half = int(len1/2)
half2 = int(len2/2)

w = WCS(hd[1].header)
w = w.dropaxis(2)
## IF TWO PARTS


w1 = w.slice(slice(0,half),slice(0,int(len2)))
w2 = w.slice(slice(half,len1),slice(0,int(len2)))
print(w1,w2)

for i in range(13):
    sli = 'a/{}a-{}_slinan.fits'.format(galaxy,i)
    hd = fits.open(sli)
    hd[1].header.update(w1.to_header())
    hd.writeto(sli,overwrite=True)
    sli = 'b/{}b-{}_slinan.fits'.format(galaxy,i)
    hdb = fits.open(sli)
    hdb[1].header.update(w2.to_header())
    hdb.writeto(sli,overwrite=True)
    print(i)


# ## IF FOUR PARTS

# w1 = w[0:half,0:half2 ]
# w4 = w[half:len1,half2:len2]

# w3 = w[0:half,half2:len2]
# w2 = w[half:len1,0:half2]
# print(w1,w2,w3,w4)

# for i in range(13):
#     sli = '{}a-{}_slinan.fits'.format(galaxy,i)
#     hd = fits.open(sli)
#     hd[1].header.update(w1.to_header())
#     hd.writeto(sli,overwrite=True)
#     sli = '{}b-{}_slinan.fits'.format(galaxy,i)
#     hdb = fits.open(sli)
#     hdb[1].header.update(w2.to_header())
#     hdb.writeto(sli,overwrite=True)
#     sli = '{}c-{}_slinan.fits'.format(galaxy,i)
#     hdc = fits.open(sli)
#     hdc[1].header.update(w3.to_header())
#     hdc.writeto(sli,overwrite=True)
#     sli = '{}d-{}_slinan.fits'.format(galaxy,i)
#     hdd = fits.open(sli)
#     hdd[1].header.update(w4.to_header())
#     hdd.writeto(sli,overwrite=True)




## IF THREE PARTS

# third = int(len1/3)

# w1 = w[0:third,:]
# w2 = w[third:third+third,:]
# w3 = w[third+third:third+third+third,:]


# print(w1,w2,w3)



# for i in range(13):
#     sli = '{}a-{}_slinan.fits'.format(galaxy,i)
#     hd = fits.open(sli)
#     hd[1].header.update(w1.to_header())
#     hd.writeto(sli,overwrite=True)
#     sli = '{}b-{}_slinan.fits'.format(galaxy,i)
#     hdb = fits.open(sli)
#     hdb[1].header.update(w2.to_header())
#     hdb.writeto(sli,overwrite=True)
#     sli = '{}c-{}_slinan.fits'.format(galaxy,i)
#     hdc = fits.open(sli)
#     hdc[1].header.update(w3.to_header())
#     hdc.writeto(sli,overwrite=True)