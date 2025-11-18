import argparse
import subprocess
import os
import sys

# --- 1. Define and Parse Arguments ---
parser = argparse.ArgumentParser(
    description="Automated pipeline for PNe detection and flux calculation using MUSE data."
)

parser.add_argument("-g", "--Galaxy", required=True, help="Galaxy name (e.g., ngc4374)")
parser.add_argument("-p", "--Path", required=True, help="Path to the main data directory (e.g., ../ngc4374/)")
parser.add_argument("-z", "--Redshift", required=True, type=float, help="Redshift of the galaxy (e.g., 0.003392)")
parser.add_argument("-sx", "--SxPath", required=True, help="Path where SExtractor default files are (e.g., /mnt/c/Users/ennis/Documents/Work/2023/pne/)")

args = parser.parse_args()

# Define common variables
GALAXY = args.Galaxy
PATH = os.path.abspath(args.Path) # Use absolute path for robustness
REDSHIFT = str(args.Redshift)
SX_PATH = args.SxPath
# Derived filenames
EMISSION_CUBE_PREFIX = f"{GALAXY}emi"
PNE_CATALOG = f"{GALAXY}_pne.cat"

# --- 2. Execution Function ---
def run_script(script_name, arguments):
    """Executes a script using subprocess and checks for errors."""
    print(f"\n--- Starting: {script_name} ---")
    command = [sys.executable, script_name] + arguments
    
    try:
        # Check=True will raise a CalledProcessError if the script returns a non-zero exit code
        subprocess.run(command, check=True, cwd=os.getcwd(), text=True, capture_output=False)
        print(f"--- Finished: {script_name} successfully ---")
    except subprocess.CalledProcessError as e:
        print(f"ERROR: {script_name} failed.")
        print(f"Return Code: {e.returncode}")
        sys.exit(1)
    except FileNotFoundError:
        print(f"ERROR: Python interpreter or script '{script_name}' not found.")
        sys.exit(1)

# --- 3. Pipeline Steps (Sequential Execution) ---

print(f"Starting PNe Pipeline for Galaxy: {GALAXY} at Path: {PATH}")

# 1. Create original continuum cube and subtract it (cont_jac.py)

original_cwd = os.getcwd()
os.chdir(PATH)
run_script('cont_jac.py', ['-g', GALAXY, '-p', '.', '-z', REDSHIFT])
os.chdir(original_cwd)
# After this, the emission lines cube is assumed to be named '{GALAXY}emi.fits'

# ---
# 2. Generate the slices from the emission lines cube (slicegenerator.py)
# This script uses the emission cube implicit to the PATH.
run_script('slicegenerator.py', ['-g', GALAXY, '-p', PATH, '-z', REDSHIFT, '-e', EMISSION_CUBE_PREFIX])

# ---
# 3. Run sextractor on the slices (sxtractor.py)
# The -p path here should point to where the slices are located (PATH).
run_script('sxtractor.py', ['-g', GALAXY, '-p', PATH, '-sx', SX_PATH])

# ---
# 4. Apply DELF method, Gaussian fitting, and flux integration (detect.py)
run_script('detect.py', ['-g', GALAXY, '-p', PATH, '-z', REDSHIFT])

# ---
# 5. Calculate final fluxes using an adaptation of Felipe Gran's method (fgran_fluxes.py)
# The -c argument uses the output catalog from detect.py.
run_script('fgran_fluxes.py', ['-g', GALAXY, '-p', PATH, '-c', PNE_CATALOG])

# ---
# 6. Generate a file with CRVAL and other values (important_data.py)
run_script('important_data.py', ['-g', GALAXY, '-p', PATH])

print("\n\n✅ Pipeline completed successfully!")