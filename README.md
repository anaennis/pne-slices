# pne-slices
General guide:

- first, masks.py creates masks, pixels have to be manually specified, useful for ppxf and can also be added to the cube for the next steps. 
- second, create_cont.py creates original continuum cube and subtracts it from the original cube, generating an emission lines cube.
- third, slicegenerator.py generates the slices from the emission lines cube
 python slicegenerator.py -g ngc4374 -p ../../ngc4374/ -z 0.003392 -e ngc4374emi
- fourth, nans.py eliminates nan values so sextractor doesn't get confused
 python nans.py -g ngc4374 -p ../../ngc4374/
- fifth, sxtractor.py runs sextractor on the slices
python sxtractor.py -g ngc4374 -p /mnt/c/Users/ennis/Documents/Work/2023/pne/ngc4374/ -sx /mnt/c/Users/ennis/Documents/Work/2023/pne/
- sixth, detect.py applies delf method
 python detect.py -g ngc4374 -p ../../ngc4374/ -z 0.003392 

- seventh, fgran_fluxes uses an adaptation of felipe gran's method to calculate the fluxes
 python fgran_fluxes.py -g ngc4374 -p ../../ngc4374/ -c ngc4374_pne.cat


- eighth, diagnostic.ipynb is a notebook that creates a match between the catalogues, allows the visualization of the lines and creates diagnostic diagrams


- extra: important_data generates a file with CRVAL and other values for the original and the emission lines cubes