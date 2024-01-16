from mpdaf.obj import Cube
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("-g","--Galaxy",help="Galaxy Name")
parser.add_argument("-p","--Path",help="Path")

args= parser.parse_args()

galaxy = args.Galaxy
path = args.Path

os.chdir(path)
nm = '{}.fits'.format(galaxy)

cube=Cube(nm)

len1 = cube.data_header['NAXIS2']
len2 = cube.data_header['NAXIS1']

half = int(len1/2)

cube1 = cube[:,0:half,:]
cube2 = cube[:,half:int(len1),:]

cube1.write('ngc4472a.fits')
cube2.write('ngc4472b.fits')
