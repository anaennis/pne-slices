import numpy as np
from astropy.io import fits
from mpdaf.obj import Spectrum,Cube,WCS,WaveCoord,iter_spe
import astropy.units as u
from numpy.polynomial import Polynomial as poly
import matplotlib.pyplot as plt


cube = Cube('ngc584.fits')

cont1 = cube.clone(data_init=np.empty, var_init=np.zeros)

count = 0
for sp, co in zip(iter_spe(cube), iter_spe(cont1)):
    try:
        co[:] = sp.poly_spec(5)
    except:
        co[:] = sp
        count += 1
print(count)

cont1.write("ngc584cont.fits")
