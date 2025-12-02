import os
import subprocess
from astropy.io import fits
import fnmatch

import argparse 

parser = argparse.ArgumentParser()

parser=argparse.ArgumentParser()

parser.add_argument("-g", "--Galaxy", help = "Galaxy name")
parser.add_argument("-p", "--Path", help = "Path where slices are")
# parser.add_argument("-sx", "--Sx", help = "Path where SExtractor default files are")


args = parser.parse_args()

galaxy = args.Galaxy
path = os.path.abspath(args.Path)

# Directory containing this script; run SExtractor commands from here as requested
THIS_DIR = os.path.dirname(os.path.abspath(__file__))

# Count available 'slinan' files in the data path
count = len(fnmatch.filter(os.listdir(path), '*nan*'))
ph = os.path.join(path, galaxy)

for i in range(count):
    i2 = str(i)
    in_file = f"{ph}-{i2}_slinan.fits[1]"
    catalog_name = f"{ph}-{i2}.cat"
    checkimg = f"{ph}check_{i2}.fits,{ph}checkrms_{i2}.fits"
    weight_img = f"{ph}-{i2}_slinan.fits[2]"

    cmd = (
        f"sex {in_file} -CATALOG_NAME {catalog_name} -SEEING_FWHM 0.8"
        f" -CHECKIMAGE_NAME {checkimg} -VERBOSE_TYPE FULL -DETECT_THRESH 1.1"
        f" -DETECT_MINAREA 4 -WEIGHT_TYPE map_Var -WEIGHT_IMAGE {weight_img}"
    )

    try:
        print(f"Running SExtractor for slice {i} (cwd={THIS_DIR}):\n  {cmd}")
        subprocess.run(cmd, shell=True, check=True, cwd=THIS_DIR, universal_newlines=True)
    except subprocess.CalledProcessError as e:
        print(f"SExtractor failed for slice {i} with return code {e.returncode}")
        # continue to next slice rather than aborting entire loop
        continue
