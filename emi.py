# code written only for when creating continuum cube with code pre january 3rd 2023 where emission line cube isnt automatically created
import numpy as np
from astropy.io import fits
from mpdaf.obj import Spectrum,Cube,WCS,WaveCoord,iter_spe
import astropy.units as u
from numpy.polynomial import Polynomial as poly
import matplotlib.pyplot as plt
import os

import argparse
import timeit

tic=timeit.default_timer()

parser=argparse.ArgumentParser()

parser.add_argument("-g", "--Galaxy", help = "Galaxy name")
parser.add_argument("-r", "--Run", help = "Folder tree")

args = parser.parse_args()

galaxy = args.Galaxy
folder = args.Run


os.chdir(folder)

nm = '{}.fits'.format(galaxy)
nm2 = '{}cont.fits'.format(galaxy)


cube = Cube(nm)

cont1 = Cube(nm2)


emi = cube - cont1
nm3 = '{}emi.fits'.format(galaxy)
emi.write(nm3)

toc = timeit.default_timer()
lap = toc - tic
print(lap)


