The pipeline consists of six sequential scripts designed to transform a raw MUSE data cube ({GALAXY}.fits) into a final PNe catalog ({GALAXY}_pne.cat).

!! Prerequisites
Before running the pipeline, ensure you have the necessary Python libraries installed (e.g., numpy, astropy, mpdaf, pandas, photutils, scipy).

You must also have the external program SExtractor (Source Extractor) installed and accessible, with a directory containing its default configuration files (e.g., .param, .conv, .mask).

If you want to run the scripts individually, you should follow this order:

1 - continuum.py 

It creates a continuum image ({GALAXY}_cont.fits).
EXAMPLE: python continuum.py -g ngc4374 -p ../ngc4374 -z 0.003392

2 - slicegenerator.py

Extracts narrow-band image slices covering the redshifted [OIII] line from the original cube.
EXAMPLE: python slicegenerator.py -g ngc4374 -p ../../ngc4374/ -z 0.003392 -e ngc4374emi

3 - emission.py
Subtracts the continuum image from all slices.
EXAMPLE: python emission.py -g ngc4374 -p ../../ngc4374 -x_start=100 y_start=100 x_end=200 y_end=200

4 - sxtractor.py

Runs SExtractor con the slices to generate preliminary source catalogs.
EXAMPLE: python sxtractor.py -g ngc4374 -p /path/to/data/ngc4374/ -sx /path/to/sextractor/defaults/

5 - detect.py

Applies the DELF method and performs Gaussian fitting on the extracted spectral trace to confirm sources and calculate total flux.
EXAMPLE: python detect.py -g ngc4374 -p ../../ngc4374/ -z 0.003392

6 - fgran_fluxes.py

Calculates aperture-corrected fluxes for the detected PNe across all slices using the generated PNe catalog.
EXAMPLE: python fgran_fluxes.py -g ngc4374 -p ../../ngc4374/ -c ngc4374_pne.cat

7 - important_data.py 

Extracts key header infromation form the original and emission cubes for later use.
EXAMPLE python important_data.py -g ngc4374 -p ../../ngc4374/


Recommended Usage: Automated Driver Script
For a cleaner and more professional execution, it is highly recommended to use the driver script, run_pipeline.py, which executes all the above steps in order:

EXAMPLE USAGE
python run_pipeline.py \
    --Galaxy ngc4374 \
    --Path /path/to/your/ngc4374/data/ \
    --Redshift 0.003392 \
    --SxPath /path/to/your/sextractor/default/files/