# pne-slices

A pipeline for detecting **planetary nebulae (PNe)** in **MUSE** integral-field data cubes and measuring their **[O III] λ5007** fluxes.

The pipeline takes a reduced MUSE cube, `{GALAXY}.fits`, and does the following:

1. Builds narrow-band "slices" around the redshifted [O III] line.
2. Subtracts the galaxy continuum from each slice.
3. Detects point sources in each slice with SExtractor.
4. Keeps only sources that appear in three consecutive slices and whose aperture spectrum is fit well by a Gaussian emission line.

The main product is a PNe catalog, `{GALAXY}_pne.cat`, with positions, integrated [O III] fluxes and line centroids. A per-wavelength aperture-photometry table, `{GALAXY}_aper_phot.csv`, is also written.

> **Citation:** If you use this code, please cite **Ennis et al. 2025**, [arXiv:2511.11479](https://arxiv.org/abs/2511.11479) (doi:[10.48550/arXiv.2511.11479](https://doi.org/10.48550/arXiv.2511.11479)).

---

## Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Input data](#input-data)
- [Quick start](#quick-start)
- [Pipeline steps](#pipeline-steps)
- [Running steps individually](#running-steps-individually)
- [Output files](#output-files)
- [Tunable parameters](#tunable-parameters)
- [Extra tools](#extra-tools)
- [Repository layout](#repository-layout)
- [Troubleshooting](#troubleshooting)
- [Acknowledgements](#acknowledgements)
- [License](#license)

---

## Requirements

**Python 3.6 or newer**, with these packages:

| Package | Used for |
|---|---|
| `numpy`, `scipy` | Array maths, Gaussian fitting, integration |
| `pandas` | Reading and writing catalogs |
| `matplotlib` | Plotting (imported by several scripts) |
| `astropy` | FITS input and output, WCS, sigma-clipped statistics, PSF model |
| `photutils` | Aperture and annulus photometry |
| `mpdaf` | Reading MUSE cubes and slicing them by wavelength |

**SExtractor (Source Extractor).** This is an external program, not a Python package. `sxtractor.py` calls it as **`sex`**, so that command must be on your `PATH`. The configuration files it needs (`default.sex`, `default.param`, `default.conv`, `default.nnw`) are included in this repository.

## Installation

```bash
git clone https://github.com/anaennis/pne-slices.git
cd pne-slices
pip install -r requirements.txt
```

To install SExtractor, use one of:

```bash
conda install -c conda-forge astromatic-source-extractor   # conda
brew install sextractor                                    # macOS / Homebrew
sudo apt install source-extractor                          # Debian / Ubuntu
```

Some packages install the program under the name `source-extractor` rather than `sex`. If so, add an alias or symlink so that `sex` works:

```bash
ln -s "$(which source-extractor)" /usr/local/bin/sex
```

## Input data

Put the MUSE cube in its own working directory. **All outputs are written to that same directory.**

```
/path/to/ngc4374/
└── ngc4374.fits      # must be named exactly {GALAXY}.fits
```

The cube must follow the standard MUSE `DATACUBE_FINAL` layout:

- **Extension 1:** `DATA`, the flux cube.
- **Extension 2:** `STAT`, the variance cube.
- **Header keywords:** `CRPIX3`, `CD3_3` and `CRVAL3`, which describe the wavelength axis.

The pipeline also needs:

- The systemic **redshift** of the galaxy.
- A rectangular **reference region**, in pixel coordinates, used to scale the continuum before subtracting it (see [Tunable parameters](#tunable-parameters)). The region should be inside the field and free of bright emission-line sources.

## Quick start

`run_pipeline.py` runs every step in order:

```bash
python run_pipeline.py \
    --Galaxy ngc4374 \
    --Path /path/to/ngc4374/ \
    --Redshift 0.003392
```

You can also pass the reference region used for continuum scaling. The values below are the defaults:

```bash
python run_pipeline.py -g ngc4374 -p /path/to/ngc4374/ -z 0.003392 \
    --x_start 200 --x_end 300 --y_start 100 --y_end 200
```

| Argument | Required | Description |
|---|---|---|
| `-g`, `--Galaxy` | yes | Galaxy name. It must match the cube filename without `.fits`. |
| `-p`, `--Path` | yes | Directory containing the cube. Outputs are written here. |
| `-z`, `--Redshift` | yes | Systemic redshift. |
| `--x_start`, `--x_end` | no | X pixel range of the reference region. Default `200`–`300`. |
| `--y_start`, `--y_end` | no | Y pixel range of the reference region. Default `100`–`200`. |

You can run `run_pipeline.py` from any directory. It finds the other scripts next to itself and runs each one with the data directory as its working directory. If any step fails, the run stops and reports which script failed.

## Pipeline steps

| # | Script | What it does | Main outputs |
|---|---|---|---|
| 1 | `continuum.py` | Sums the cube over the rest-frame window 5031–5154 Å, shifted by the redshift, to make a continuum image. | `{G}_cont.fits` |
| 2 | `slicegenerator.py` | Extracts wavelength slices around the redshifted [O III] λ5007 line, over a window of ±0.0016 in *z* (about ±480 km/s). It writes 16 single-channel slices and 13 slices that each sum three channels. | `{G}-{i}_sliun.fits` (i = 0–15)<br>`{G}-{i}_sli.fits` (i = 0–12) |
| 3 | `emission.py` | **(a)** Scales the continuum image to each single-channel slice, using the mean flux ratio in the reference region, and subtracts it. **(b)** Adds each run of three consecutive subtracted slices together, combining the errors in quadrature, to make the detection images. | `sub-{G}-{i}.fits` (i = 0–15)<br>`{G}-{i}_slinan.fits` (i = 0–12) |
| 4 | `sxtractor.py` | Runs SExtractor on each detection image (extension 1), using the variance (extension 2) as a weight map. | `{G}-{i}.cat`<br>`{G}check_{i}.fits`, `{G}checkrms_{i}.fits` |
| 5 | `detect.py` | Implements the DELF method. It keeps sources that appear in three consecutive SExtractor catalogs within 0.8 px, removes duplicates, and measures 3 px aperture photometry on all 16 continuum-subtracted slices to build a spectrum for each source. It fits a Gaussian to each spectrum and keeps sources that pass quality cuts. The flux is the integral of the fitted Gaussian. | `{G}_pne.cat`<br>`{G}_vrad_wcs.cat` |
| 6 | `fgran_fluxes.py` | Measures 5 px aperture photometry at each PN position in **every wavelength channel of the original cube**. It subtracts a sigma-clipped annulus background (5–8 px) and applies a Moffat-PSF aperture correction. | `{G}_aper_phot.csv` |
| 7 | `important_data.py` | Saves the wavelength-axis header keywords (`CRPIX3`, `CD3_3`, `CRVAL3`) so that channel numbers can be converted to wavelengths later. | `{G}_dat.cat` |

`{G}` stands for the galaxy name.

## Running steps individually

Each script can also be run by itself, which helps when re-running or debugging one stage. Run them **in the order below**, because each step reads the files written by the steps before it.

```bash
# 1. Continuum image
python continuum.py -g ngc4374 -p /path/to/ngc4374/ -z 0.003392

# 2. Narrow-band slices around [O III]
python slicegenerator.py -g ngc4374 -p /path/to/ngc4374/ -z 0.003392

# 3. Continuum subtraction + 3-slice co-addition
python emission.py -g ngc4374 -p /path/to/ngc4374/ \
    --x_start 200 --x_end 300 --y_start 100 --y_end 200

# 4. SExtractor source detection
python sxtractor.py -g ngc4374 -p /path/to/ngc4374/

# 5. DELF matching + Gaussian line fitting → PNe catalog
python detect.py -g ngc4374 -p /path/to/ngc4374/ -z 0.003392

# 6. Aperture-corrected photometry across the full cube
python fgran_fluxes.py -g ngc4374 -p /path/to/ngc4374/ -c ngc4374_pne.cat

# 7. Wavelength-axis header information
python important_data.py -g ngc4374 -p /path/to/ngc4374/
```

Notes:

- Each script changes into the `-p` directory before it reads or writes anything. Relative file arguments, such as `-c` in step 6, are relative to that directory.
- `sxtractor.py` runs SExtractor from the **repository directory**, so that it picks up the `default.*` configuration files there. Data paths are passed to SExtractor as absolute paths, so the location you run it from does not matter.

## Output files

### `{GALAXY}_pne.cat`: the final PNe catalog

A whitespace-separated table with no header row. There is one row per PN, sorted by `x` and then `y`.

| Column | Description |
|---|---|
| `ra` | Right ascension, J2000, in degrees (from SExtractor `ALPHA_J2000`) |
| `dec` | Declination, J2000, in degrees (from SExtractor `DELTA_J2000`) |
| `x`, `y` | Pixel position in the cube |
| `flux` | Integrated [O III] flux: the fitted Gaussian integrated over the slice window, in cube flux units × Å |
| `cen` | Fitted line centre, in Å |

To read it:

```python
import pandas as pd
pne = pd.read_csv("ngc4374_pne.cat", sep=r"\s+", names=["ra", "dec", "x", "y", "flux", "cen"])
```

### Other products

| File | Contents |
|---|---|
| `{G}_vrad_wcs.cat` | `ra dec vrad`, with `vrad` in km/s computed from the fitted line centre. It has one row for **every three-slice candidate**, including those rejected by the fit cuts. Rejected candidates have `vrad = 0`. |
| `{G}_aper_phot.csv` | One row per PN, with columns `id`, `xcenter`, `ycenter`, then `aper_sum_bkgsub_N` and `total_error_N` for each channel `N` of the original cube, then the aperture-corrected fluxes `aper_sum_bkgsub_N_corrected`. |
| `{G}_dat.cat` | Eight lines: the cube filename, `CRPIX3`, `CD3_3`, `CRVAL3`, then the same set labelled with the `{G}emi.fits` filename. Both sets are read from `{G}.fits`. Convert channel `N` to wavelength with λ = (N + 1 − CRPIX3) · CD3_3 + CRVAL3. |
| `{G}_cont.fits` | Continuum image |
| `{G}-{i}_sliun.fits` | Single-channel slices |
| `{G}-{i}_sli.fits` | Three-channel slices taken straight from the cube, without continuum subtraction. Later steps do not use these. |
| `sub-{G}-{i}.fits` | Continuum-subtracted single-channel slices. HDU 1 holds the data and HDU 2 the variance. |
| `{G}-{i}_slinan.fits` | Continuum-subtracted three-channel detection images. HDU 1 holds the data and HDU 2 the combined error. |
| `{G}-{i}.cat` | SExtractor catalog for each detection image. Its columns are set by `default.param`. |
| `{G}check_{i}.fits`, `{G}checkrms_{i}.fits` | SExtractor background and background-RMS check images |

## Tunable parameters

Most settings are fixed in the scripts and chosen for MUSE Wide Field Mode data (0.2″ pixels, 1.25 Å channels). Change them in the files listed below if your data differ.

| Setting | Value | Where |
|---|---|---|
| Continuum window, rest frame | 5031–5154 Å | `continuum.py` |
| [O III] reference wavelength | 5008 Å | `slicegenerator.py`, `detect.py` |
| Velocity window | *z* ± 0.0016 (about ±480 km/s) | `slicegenerator.py`, `detect.py` |
| Number of single-channel / three-channel slices | 16 / 13 | `slicegenerator.py`, `emission.py` |
| Continuum-scaling reference region | x 200–300, y 100–200 px | CLI flags of `emission.py` / `run_pipeline.py` |
| SExtractor overrides | `DETECT_THRESH 1.1`, `DETECT_MINAREA 4`, `SEEING_FWHM 0.8`, `WEIGHT_TYPE MAP_VAR` | `sxtractor.py` |
| Other SExtractor settings | `PIXEL_SCALE 0.2`, `MAG_ZEROPOINT 32.2`, `BACK_SIZE 64`, … | `default.sex` |
| Distance tolerance for matching across slices | 0.8 px | `detect.py` |
| Aperture radius for the spectrum | 3 px | `detect.py` |
| Channel width | 1.25 Å | `detect.py` |
| Gaussian-fit quality cuts | amplitude S/N ≥ 0.5, reduced χ² ≤ 2, 0.001 < FWHM < 10 Å | `detect.py` |
| Radius for removing duplicates | 0.5 px | `detect.py` |
| Final aperture / background annulus | r = 5 px / 5–8 px | `fgran_fluxes.py` |
| PSF for aperture correction | Moffat, γ = 3, α = 2.45 | `fgran_fluxes.py` |

`continuum.py` and `slicegenerator.py` read the wavelength axis from spaxel (100, 100). The default reference region also assumes the cube is at least 300 × 200 spaxels. For smaller cubes, choose a reference region that fits inside the field.

## Extra tools

These are **not** run by `run_pipeline.py`. They are included as optional helpers.

- **`crop.py`**: splits a cube into two halves along the y axis, which helps when a cube is too large to process in one go. The output names are currently fixed as `ngc4472a.fits` and `ngc4472b.fits`. Edit them in the script before using it on another galaxy.

  ```bash
  python crop.py -g ngc4472 -p /path/to/ngc4472/
  ```

- **`fix_wcs.py`**: after each half of a split cube has been run through the pipeline separately, in sub-folders `a/` and `b/`, this script writes the correct WCS from a full-field reference file (`stitched.fits`) into each half's `_slinan.fits` detection images. Commented-out blocks in the script handle splits into three or four parts.

  ```bash
  python fix_wcs.py -g ngc4472 -p /path/to/slices/
  ```

- **`diagnostic-NEW-ngc1387.ipynb`**: a Jupyter notebook, provided for reference, with the diagnostic checks and plots made for NGC 1387. It is not needed to run the pipeline.

## Repository layout

```
pne-slices/
├── run_pipeline.py        # Runs every step in order
├── continuum.py           # Step 1: continuum image
├── slicegenerator.py      # Step 2: [O III] narrow-band slices
├── emission.py            # Step 3: continuum subtraction + co-addition
├── sxtractor.py           # Step 4: SExtractor detection
├── detect.py              # Step 5: DELF matching + Gaussian fitting
├── fgran_fluxes.py        # Step 6: aperture-corrected photometry
├── important_data.py      # Step 7: wavelength-axis header info
├── crop.py                # Extra: split a large cube in two
├── fix_wcs.py             # Extra: restore WCS on split-cube slices
├── default.sex            # SExtractor configuration
├── default.param          # SExtractor output columns
├── default.conv           # SExtractor detection filter
├── default.nnw            # SExtractor star/galaxy neural network
├── diagnostic-NEW-ngc1387.ipynb
├── requirements.txt
└── LICENSE
```

## Troubleshooting

- **`sex: command not found`**: SExtractor is either not installed or installed as `source-extractor`. See [Installation](#installation).
- **SExtractor runs on the wrong files, or on too many files.** `sxtractor.py` counts the files in the data directory whose names contain `nan` to work out how many detection images there are. Keep other files with `nan` in their names out of that directory.
- **`detect.py` writes an empty catalog.** No source appeared in three consecutive slices. Check the redshift. Also check that the reference region in `emission.py` is free of emission and gives a sensible continuum subtraction: the scale factor `k` is printed for each slice.
- **`fgran_fluxes.py` is slow.** It measures photometry in every channel of the full cube, about 3,700 channels for MUSE, so it takes a while for large catalogs.
- **Re-running.** Outputs are written to the data directory and are overwritten on each run.

## Acknowledgements

The photometry in `fgran_fluxes.py` is adapted from a method by Felipe Gran, which builds on the [photutils](https://photutils.readthedocs.io/en/stable/aperture.html) aperture-photometry routines.

If you use this pipeline in published work, please cite **Ennis et al. 2025** ([arXiv:2511.11479](https://arxiv.org/abs/2511.11479)).

## License

This project is released under the [MIT License](LICENSE).
