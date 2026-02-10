"""
create_emission_images.py: Background Subtraction and Co-addition

This script performs two main functions:
1. Background Subtraction: Subtracts a scaled continuum image from each unbinned 
   emission-line slice (e.g., '{GALAXY}-{i}_sliun.fits'). The scaling factor 
   (k) is determined by the mean flux ratio in a user-defined sky region.
2. Co-addition: Combines three consecutive background-subtracted slices to 
   create the final binned emission-line images (e.g., '{GALAXY}-{i}_slinan.fits'),
   propagating errors in quadrature.

This process prepares the final images for SExtractor detection.

Requires: astropy, numpy, os, argparse
"""


import argparse
import os
import numpy as np
from astropy.io import fits

def process_single_slice(galaxy, cont_filename, slice_number, region_slice):
    """
    Processes a single unbinned slice by performing background subtraction 
    using a scaled continuum image.

    :param galaxy: Name of the galaxy.
    :type galaxy: str
    :param cont_filename: The filename of the continuum image.
    :type cont_filename: str
    :param slice_number: The index of the unbinned slice to process (0-based).
    :type slice_number: int
    :param region_slice: NumPy slice object (np.s_) defining the background region.
    :type region_slice: numpy.s_
    :returns: True if successful, False otherwise.
    :rtype: bool
    """
    emission_line_filename = f'{galaxy}-{slice_number}_sliun.fits'
    output_filename = f'sub-{galaxy}-{slice_number}.fits'

    try:
        # Open FITS files for the unbinned slice and the continuum image
        with fits.open(emission_line_filename) as hdus, \
             fits.open(cont_filename) as hduc:

            emission_data = hdus[1].data
            
            cont_data = hduc[1].data / 101 
            
            # --- Calculate Scaling Factor (k) ---
            # Mean flux in the defined background region for continuum image
            ctd = np.nanmean(cont_data[region_slice])
            # Mean flux in the defined background region for emission slice
            std = np.nanmean(emission_data[region_slice])

            # Calculate k: the scaling factor
            if ctd == 0:
                print(f"Warning: Continuum mean is zero in region for slice {slice_number}. Setting k=1.")
                k = 1.0
            else:
                k = std / ctd
            
            # --- Perform Background Subtraction ---
            resta = emission_data - k * cont_data

            # Update the data in the HDU list
            hdus[1].data = resta.astype(np.float32)
            
            # Copy WCS from original slice to output headers
            # Copy WCS keywords from the original image HDU header to all output HDUs
            if hdus[1].header:
                for hdu_idx in [0, 1, 2]:  # primary, data, error
                    if hdu_idx < len(hdus):
                        # Copy WCS-related keywords
                        wcs_keywords = ['WCSAXES', 'CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2', 
                                       'CTYPE1', 'CTYPE2', 'CUNIT1', 'CUNIT2',
                                       'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
                                       'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2',
                                       'CDELT1', 'CDELT2', 'CROTA1', 'CROTA2']
                        for key in wcs_keywords:
                            if key in hdus[1].header:
                                hdus[hdu_idx].header[key] = hdus[1].header[key]
            
            # --- Write the result ---
            hdus.writeto(output_filename, overwrite=True)
            print(f"✅ Subtracted: {emission_line_filename} (k={k:.4f}) -> {output_filename}")
            return True

    except FileNotFoundError:
        print(f"❌ Error: Required file not found for slice {slice_number}. Skipping.")
        return False
    except Exception as e:
        print(f"An error occurred while processing slice {slice_number}: {e}")
        return False


def generate_coadded_images(galaxy, total_slices):
    """
    Co-adds three consecutive background-subtracted slices, propagating errors 
    in quadrature, and saves the final binned images.

    :param galaxy: Name of the galaxy.
    :type galaxy: str
    :param total_slices: The total number of unbinned 'sub-*.fits' files available.
    :type total_slices: int
    :returns: None
    """
    print("\n--- Starting Co-addition (3-slice binning) ---")

    # Dictionaries to store loaded FITS data and error data
    fits_data = {}
    fits_errors = {}
    
    # 1. Load all background-subtracted slices
    for i in range(total_slices):
        filename = f"sub-{galaxy}-{i}.fits"
        try:
            with fits.open(filename) as hdul:
                if len(hdul) >= 3:

                    fits_data[i] = hdul[1].data
                    fits_errors[i] = hdul[2].data
                    if i == 0:
                         print(f"Loaded background-subtracted files (expecting prefix 'sub-{galaxy}-')")
                else:
                    print(f"Warning: {filename} missing HDU 2 (Error). Skipping slice {i}.")
        except FileNotFoundError:
            continue 

    # 2. Generate 3-slice coadded images

    num_coadds = 13
    
    try:
        header_template_filename = f'{galaxy}-1_sliun.fits'
        with fits.open(header_template_filename) as h:
             # Use the primary header for the new file structure
             template_header = h[0].header.copy() 
    except FileNotFoundError:
        print(f"Fatal Error: Could not load header template from {header_template_filename}. Cannot create final FITS files.")
        return

    for i in range(1, num_coadds + 1):
        # Coadd i is centered at index (i+1) and uses slices (i, i+1, i+2)
        slice_indices = [i, i + 1, i + 2] 
        
        # Check if all 3 required slices are available
        if all(idx in fits_data for idx in slice_indices):
            
            # 3. Coadd Data (Sum)
            data_to_coadd = [fits_data[idx] for idx in slice_indices]
            coadded_data = np.sum(data_to_coadd, axis=0)

            # 4. Propagate Errors (Quadrature Sum for Sum)
            error_squares = sum(fits_errors[idx]**2 for idx in slice_indices)
            coadded_error = np.sqrt(error_squares)

            coadded_data = coadded_data.astype(np.float32)
            coadded_error = coadded_error.astype(np.float32)

            # 5. Save the Coadded Image
            # Output format: HDU 1 = Data, HDU 2 = Error
            
            output_index = i - 1
            output_filename = f"{galaxy}-{output_index}_slinan.fits"
            
            # Load WCS from one of the center slice's original file
            center_slice_filename = f'{galaxy}-{slice_indices[1]}_sliun.fits'
            wcs_header = None
            try:
                with fits.open(center_slice_filename) as wcs_hdul:
                    if len(wcs_hdul) > 1 and wcs_hdul[1].header:
                        wcs_header = wcs_hdul[1].header.copy()
            except:
                pass
            
            # Create a new FITS HDUList
            primary_hdu = fits.PrimaryHDU(header=template_header) 
            image_hdu = fits.ImageHDU(coadded_data, name='DATA')
            error_hdu = fits.ImageHDU(coadded_error, name='ERR')
            
            # Copy WCS keywords to all HDUs if available
            if wcs_header:
                wcs_keywords = ['WCSAXES', 'CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2', 
                               'CTYPE1', 'CTYPE2', 'CUNIT1', 'CUNIT2',
                               'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
                               'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2',
                               'CDELT1', 'CDELT2', 'CROTA1', 'CROTA2']
                for hdu in [primary_hdu, image_hdu, error_hdu]:
                    for key in wcs_keywords:
                        if key in wcs_header:
                            hdu.header[key] = wcs_header[key]
            
            hdul = fits.HDUList([primary_hdu, image_hdu, error_hdu])
            
            hdul.writeto(output_filename, overwrite=True)
            print(f"✅ Coadded: Slices {slice_indices} -> {output_filename}")
            hdul.close()

        else:
            print(f"Warning: Missing slices for coadd centered at index {slice_indices[1]}. Skipping coadd {i}.")


def main():
    """Main execution function to run the two-part process."""
    parser = argparse.ArgumentParser(
        description="Creates final emission-line images by performing continuum subtraction and co-addition."
    )
    parser.add_argument("-g", "--Galaxy", required=True, help="Galaxy name (e.g., ngc1387).")
    parser.add_argument("-p", "--Path", required=True, help="Path to the data directory.")
    parser.add_argument("--x_start", type=int, default=200, help="X-coordinate start for background region.")
    parser.add_argument("--x_end", type=int, default=300, help="X-coordinate end for background region.")
    parser.add_argument("--y_start", type=int, default=100, help="Y-coordinate start for background region.")
    parser.add_argument("--y_end", type=int, default=200, help="Y-coordinate end for background region.")
    
    args = parser.parse_args()

    # --- Setup ---
    GALAXY = args.Galaxy
    CONT_FILENAME = '{}_cont.fits'.format(GALAXY)
    TOTAL_SLICES = 16
    
    # Define the background region based on input coordinates
    REGION_SLICE = np.s_[args.y_start:args.y_end, args.x_start:args.x_end]
    
    if not os.path.exists(args.Path):
        print(f"Error: Path not found: {args.Path}")
        return
    os.chdir(args.Path)
    
    # --- Part 1: Background Subtraction ---
    print("--- Part 1: Processing individual slices for background subtraction ---")
    
    # We iterate from 0 up to (but not including) TOTAL_SLICES
    for i in range(TOTAL_SLICES):
        process_single_slice(GALAXY, CONT_FILENAME, i, REGION_SLICE)

    # --- Part 2: Co-addition and Error Propagation ---
    generate_coadded_images(GALAXY, TOTAL_SLICES)
    
    print("\n✅ Emission line images created successfully.")

if __name__ == "__main__":
    main()