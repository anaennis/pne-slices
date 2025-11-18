from mpdaf.obj import Cube,Spectrum
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--Galaxy", help="Galaxy name")
parser.add_argument("-p", "--Path", help="Path")
parser.add_argument("-z", "--Redshift", help="Redshift")

args = parser.parse_args()

galaxy = args.Galaxy
path = args.Path
z = args.Redshift
z = float(z)
os.chdir(path)
nm = '{}.fits'.format(galaxy)

cb = Cube(nm)
sp = cb[:,100,100]

oiii_1 = 5031.0*(1.0+z)
oiii_2 = 5154.0*(1.0+z)
ka,kb = sp.wave.pixel([oiii_1,oiii_2],nearest=True)

tmp = cb[ka:kb,:,:]

im = tmp.sum(axis=0)
nm_cont = '{}_cont.fits'.format(galaxy)
im.write(nm_cont)
