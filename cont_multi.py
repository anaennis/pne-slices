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


cube = Cube(nm)



cont2 = cube.loop_spe_multiprocessing(f=Spectrum.poly_specb,cpu=48,deg=5)

nm2 = '{}cont.fits'.format(galaxy)

cont2.write(nm2)

toc = timeit.default_timer()
lap = toc - tic
print(lap)