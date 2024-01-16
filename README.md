# pne-slices
General guide:

- first, masks.py creates masks, pixels have to be manually specified, useful for ppxf and can also be added to the cube for the next steps. 
- second, create_cont.py creates original continuum cube and subtracts it from the original cube, generating an emission lines cube.
- third, slicegenerator.py generates the slices from the emission lines cube
- fourth, nans.py eliminates nan values so sextractor doesn't get confused
- fifth, sxtractor.py runs sextractor on the slices
- sixth, detect.py applies delf method
- seventh, fgran_fluxes uses an adaptation of felipe gran's method to calculate the fluxes
- eighth, diagnostic.ipynb is a notebook that creates a match between the catalogues, allows the visualization of the lines and creates diagnostic diagrams
- extra: important_data generates a file with CRVAL and other values for the original and the emission lines cubes
