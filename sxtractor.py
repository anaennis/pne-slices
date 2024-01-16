import os
from astropy.io import fits
import fnmatch

import argparse 

parser = argparse.ArgumentParser()

parser=argparse.ArgumentParser()

parser.add_argument("-g", "--Galaxy", help = "Galaxy name")
parser.add_argument("-p", "--Path", help = "Path where slices are")
parser.add_argument("-sx", "--Sx", help = "Path where SExtractor default files are")


args = parser.parse_args()

galaxy = args.Galaxy
path = args.Path
path2 = args.Sx

os.chdir(path)
count = len(fnmatch.filter(os.listdir(path), '*nan*'))
ph = path + "/" + galaxy

os.chdir(path2)

i = 0
for i in range(count):
    i2 = str(i)
    pt1 = "sextractor "+ph+"-"+i2+"_slinan.fits[1]"+" -CATALOG_NAME "+ph+"-"+i2+".cat -SEEING_FWHM 0.8 -CHECKIMAGE_NAME "+ph+"check_"+i2+".fits,"+ph+"checkrms_"+i2+".fits"
    pt2 = " -VERBOSE_TYPE FULL -DETECT_THRESH 1.1 -DETECT_MINAREA 4 -WEIGHT_TYPE map_Var -WEIGHT_IMAGE "+ph+"-"+i2+"_slinan.fits[2]"
    sx = pt1+pt2
    os.system(sx)