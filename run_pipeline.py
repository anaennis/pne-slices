import argparse
import subprocess
import os
import sys
import glob

# --- 1. Define and Parse Arguments ---
parser = argparse.ArgumentParser(
    description="Automated pipeline for PNe detection and flux calculation using MUSE data."
)

parser.add_argument("-g", "--Galaxy", required=True, help="Galaxy name (e.g., ngc4374)")
parser.add_argument("-p", "--Path", required=True, help="Path to the main data directory (e.g., ../ngc4374/)")
parser.add_argument("-z", "--Redshift", required=True, type=float, help="Redshift of the galaxy (e.g., 0.003392)")
# Optional region arguments to forward to emission.py (defaults match emission.py)
parser.add_argument("--x_start", type=int, default=200, help="X-coordinate start for background region.")
parser.add_argument("--x_end", type=int, default=300, help="X-coordinate end for background region.")
parser.add_argument("--y_start", type=int, default=100, help="Y-coordinate start for background region.")
parser.add_argument("--y_end", type=int, default=200, help="Y-coordinate end for background region.")

args = parser.parse_args()

# Define common variables
GALAXY = args.Galaxy
PATH = os.path.abspath(args.Path) # Use absolute path for robustness
REDSHIFT = str(args.Redshift)
# Determine where helper scripts live. Scripts live in the same folder as this runner.
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = THIS_DIR
# Derived filenames
EMISSION_CUBE_PREFIX = f"{GALAXY}emi"
PNE_CATALOG = f"{GALAXY}_pne.cat"

# --- 2. Execution Function ---
def run_script(script_name, arguments):
    """Executes a script using subprocess and checks for errors."""
    print(f"\n--- Starting: {script_name} ---")
    # Resolve script path inside the SCRIPTS_DIR unless an absolute path was provided
    if os.path.isabs(script_name):
        script_path = script_name
    else:
        script_path = os.path.join(SCRIPTS_DIR, script_name)
    # If script_path doesn't exist, try to locate it anywhere under THIS_DIR
    if not os.path.isfile(script_path):
        print(f"Warning: script not found at expected path: {script_path}")
        # try a recursive search for the script name in the repository
        candidates = []
        for fn in glob.glob(os.path.join(THIS_DIR, '**', script_name), recursive=True):
            if os.path.isfile(fn):
                candidates.append(fn)
        if len(candidates) == 1:
            script_path = candidates[0]
            print(f"Found script at: {script_path}; will use this path.")
        elif len(candidates) > 1:
            print("Multiple candidates found for script; please disambiguate:")
            for c in candidates:
                print(f" - {c}")
            sys.exit(1)
        else:
            # show available files in scripts dir for easier debugging
            print(f"Looking in scripts directory: {SCRIPTS_DIR}")
            try:
                for p in sorted(os.listdir(SCRIPTS_DIR)):
                    print(f" - {p}")
            except Exception:
                pass
            print(f"ERROR: Script '{script_name}' not found. Ensure it exists next to run_pipeline.py or provide an absolute path.")
            sys.exit(1)

    if not os.path.exists(sys.executable):
        print(f"ERROR: Python interpreter not found at '{sys.executable}'.")
        sys.exit(1)

    command = [sys.executable, script_path] + arguments

    try:
        # Run the script with the working directory set to the data/output path (PATH)
        # so scripts operate on data located elsewhere.
        print(f"Running: {command}  (cwd={PATH})")
        # Use `universal_newlines=True` for Python 3.6 compatibility
        subprocess.run(command, check=True, cwd=PATH, universal_newlines=True)
        print(f"--- Finished: {script_name} successfully ---")
    except subprocess.CalledProcessError as e:
        print(f"ERROR: {script_name} failed.")
        print(f"Return Code: {e.returncode}")
        sys.exit(1)

# --- 3. Pipeline Steps (Sequential Execution) ---

print(f"Starting PNe Pipeline for Galaxy: {GALAXY} at Path: {PATH}")
print(f"Resolved scripts directory: {SCRIPTS_DIR}")

# 1. Create original continuum cube and subtract it (continuum.py)

# Run continuum.py (scripts live in `scripts/`, data in `PATH`).
# Pass the full data path to the script using `-p` so it knows where to write outputs.
run_script('continuum.py', ['-g', GALAXY, '-p', PATH, '-z', REDSHIFT])
# After this, the emission lines cube is assumed to be named '{GALAXY}emi.fits' in `PATH`

# ---
# 2. Generate the slices from the emission lines cube (slicegenerator.py)
# This script uses the emission cube implicit to the PATH.
run_script('slicegenerator.py', ['-g', GALAXY, '-p', PATH, '-z', REDSHIFT])

# ---
# 3. Generate the emission only slices (emission.py)
# The -p path here should point to where the slices are located (PATH).
run_script('emission.py', [
    '-g', GALAXY,
    '-p', PATH,
    '--x_start', str(args.x_start),
    '--x_end', str(args.x_end),
    '--y_start', str(args.y_start),
    '--y_end', str(args.y_end),
])

# ---
# 4. Run sextractor on the slices (sxtractor.py)
# The -p path here should point to where the slices are located (PATH).
run_script('sxtractor.py', ['-g', GALAXY, '-p', PATH])

# ---
# 5. Apply DELF method, Gaussian fitting, and flux integration (detect.py)
run_script('detect.py', ['-g', GALAXY, '-p', PATH, '-z', REDSHIFT])

# ---
# 6. Calculate final fluxes using an adaptation of Felipe Gran's method (fgran_fluxes.py)
# The -c argument uses the output catalog from detect.py.
run_script('fgran_fluxes.py', ['-g', GALAXY, '-p', PATH, '-c', PNE_CATALOG])

# ---
# 7. Generate a file with CRVAL and other values (important_data.py)
run_script('important_data.py', ['-g', GALAXY, '-p', PATH])

print("\n\n✅ Pipeline completed successfully!")