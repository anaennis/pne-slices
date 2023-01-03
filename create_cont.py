import numpy as np
from astropy.io import fits
from mpdaf.obj import Spectrum,Cube,WCS,WaveCoord,iter_spe
import astropy.units as u
from numpy.polynomial import Polynomial as poly
import matplotlib.pyplot as plt


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

cont1 = cube.clone(data_init=np.empty, var_init=np.zeros)

count = 0
for sp, co in zip(iter_spe(cube), iter_spe(cont1)):
    try:
        co[:] = sp.poly_spec(5)
    except:
        co[:] = sp
        count += 1
print(count)

nm2 = '{}cont.fits'.format(galaxy)

cont1.write(nm2)

toc = timeit.default_timer()
lap = toc - tic
print(lap)
