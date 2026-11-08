# INDO36 — NEMO 4.2.2 preprocessing workflow

**Created by Dr Dwiyoga Nugroho, Toulouse, November 2026**

Python tools for preparing the INDO36 regional NEMO model: ocean initial and open-boundary conditions, tidal harmonics and self-attraction/loading (SAL), ERA5 forcing and interpolation weights, river discharge, chlorophyll interpolation weights, and a generated NEMO namelist.

> **Status:** This guide is based on the 12 supplied Python files (11 numbered scripts and `config.py`). It documents their implemented interfaces, not a validated NEMO production run. The uploaded set does **not** include the AMM12 namelist template, input NetCDF data, tidal atlases, NEMO executable, or WEIGHTS executables. These must be obtained separately. Some forcing/namelist compatibility issues remain to be checked before running NEMO.

## 1. Requirements

- Linux, HPC cluster, or Windows Subsystem for Linux (WSL2); sufficient memory and disk space for the 2152 × 1724 × 75 example NEMO grid.
- Conda/Mamba, Python 3.11; internet access and account credentials for the download steps.
- A **separately built NEMO 4.2.2 executable**, compatible domain configuration and enabled physics; this Python environment does not install/compile NEMO.
- **Compiled NEMO `tools/WEIGHTS` executables** `scripgrid.exe`, `scrip.exe`, `scripshape.exe` for scripts 8 and 10. Conda does not install these project-specific Fortran executables.
- Optional command-line utilities `ncdump`, `ncks` (NetCDF/NCO) for inspecting files.

### Create the Python environment

```bash
conda env create -f environment.yml
conda activate indo36
python -c "import numpy, scipy, pandas, xarray, netCDF4, cdsapi, copernicusmarine; print('Dependencies OK')"
```

If the environment already exists:

```bash
conda env update -n indo36 -f environment.yml --prune
```

**File naming:** `environment.yml` is the recommended Conda convention. `environment.nml` in this repository is an identical YAML-formatted copy, **not** a Fortran namelist. You may also run `conda env create -f environment.nml`, because Conda reads the contents rather than requiring a `.yml` extension.

### Download account setup

- **Copernicus Marine** (`1_download_cmems.py`): authenticate using the `copernicusmarine` CLI or its supported credentials configuration; register at https://marine.copernicus.eu/.
- **Copernicus Climate Data Store** (`2_download_era5.py`): configure the `cdsapi` credentials and accept the ERA5 dataset terms at https://cds.climate.copernicus.eu/.
- **Copernicus Emergency Management / EWDS** (`3_download_glofass.py`): configure credentials for `https://ewds.climate.copernicus.eu/api` and accept the relevant GloFAS dataset terms at https://ewds.climate.copernicus.eu/ . The code explicitly uses the EWDS endpoint, which may require its own credentials/configuration.

Never commit API keys or credentials to GitHub.

## 2. Repository structure

```text
INDO36/
├── README.md
├── environment.yml
├── environment.nml
├── config.py
├── 1_download_cmems.py
├── 2_download_era5.py
├── 3_download_glofass.py
├── 4_make_boundary_coordinates.py
├── 5_make_initial.py
├── 6_make_tides.py
├── 7_make_ocean_boundary.py
├── 8_make_era5_weights.py
├── 9_make_river_runoff.py
├── 10_make_chlorophyll.py
├── 11_make_namelist_cfg.py
└── namelist_cfg_AMM12_template   # REQUIRED, not supplied with uploads
```

All numbered scripts use settings from `config.py`, except script 11 dynamically loads the local `config.py` and reads `namelist_cfg_AMM12_template` relative to the current working directory. Run scripts from the repository directory.

## 3. Basic data directory and required files

The supplied `config.py` assumes:

```text
/mnt/d/indo36/
├── Program/                       # Python scripts, config.py, template
├── Data/
│   ├── coordinates.nc              # REQUIRED: NEMO horizontal coordinates
│   ├── mesh_mask.nc                # REQUIRED: NEMO T/U/V/F wet masks and levels
│   ├── domain_cfg.nc               # REQUIRED: NEMO domain geometry/depth/metrics
│   ├── GLORYS12V1_1dAV_20190101_20190102_gridT_R20190102.nc
│   ├── GLORYS12V1_1dAV_20190101_20190102_gridS_R20190102.nc
│   ├── GLORYS12V1_1dAV_20190101_20190102_gridU_R20190102.nc
│   ├── GLORYS12V1_1dAV_20190101_20190102_gridV_R20190102.nc
│   ├── GLORYS12V1_1dAV_20190101_20190102_grid2D_R20190102.nc
│   ├── CMEMS/                      # Downloaded CMEMS daily five-grid files
│   ├── TPXO10/
│   │   ├── h_tpxo10.v2a.nc         # TPXO10 sea-level harmonic atlas
│   │   ├── u_tpxo10.v2a.nc         # TPXO10 transport harmonic atlas
│   │   ├── grid_tpxo10v2a.nc      # TPXO10 bathymetry/masks
│   │   └── sal_tpxo10.v2a.nc      # REQUIRED if TIDE_GENERATE_SAL=True
│   ├── FES2014/                    # Alternative tidal atlas; see below
│   └── merged_ESACCI_BIOMER4V1R1_CHL_REG05.nc
└── run/
    ├── log/
    ├── era5/
    ├── glofass_data/
    ├── tides/
    ├── bdy/
    ├── runoff/
    └── weights/
```

These are example paths from the uploaded configuration, not a claim that the files exist on your computer. Adjust `DATA_DIR`, `RUN_DIR`, and other paths in `config.py` for your machine or cluster. **Do not commit the large input datasets to GitHub**; document their acquisition and licensing separately.

### NEMO grid files

| File | Used by | Essential contents / checks |
|---|---|---|
| `coordinates.nc` | 5, 6, 8, 10 | Curvilinear NEMO T/U/V coordinates (`glamt/gphit`, `glamu/gphiu`, `glamv/gphiv`), matching model horizontal dimensions |
| `mesh_mask.nc` | 4, 5, 7, 9 | `tmask`, `umask`, `vmask`, `fmask`; vertical level information; script 5 detects the level count automatically |
| `domain_cfg.nc` | 5, 7, 9 | Vertical depths and thicknesses (`gdept_0`, `e3t_0`, etc.), horizontal metrics (`e1t`, `e2t`), and bathymetry/wet-domain fields |
| `coordinates.bdy.nc` | 6, 7 | **Generated by script 4**, using `mesh_mask.nc`; boundary positions, rim indices, and geographic coordinates |

For the example grid, the scripts report **2152 × 1724** horizontal cells and **75** levels. Grid files must be consistent with the actual compiled NEMO domain. An ORCA36 native grid crop is not by itself a complete validated `mesh_mask.nc`/`domain_cfg.nc` pair.

### Ocean sources: GLORYS versus CMEMS

- **GLORYS native**: script 5 expects one first-day file each for `gridT`, `gridS`, `gridU`, `gridV`, `grid2D` in `GLORYS_INITIAL_DIR`. Its variable aliases include `votemper`, `vosaline`, `vozocrtx`, `vomecrty`, `sossheig`; native source-grid geometry is read from the source files.
- **CMEMS**: script 1 downloads daily mean `thetao`, `so`, `uo`, `vo`, `zos`, renames to the aliases above, and saves five daily NetCDF files per day in `CMEMS_DOWNLOAD_DIR`. This branch uses geographic lon/lat source axes; script 5 rotates geographic velocity into NEMO grid axes.
- Set **`INITIAL_SOURCE_TYPE`** and **`BOUNDARY_SOURCE_TYPE`** independently to `GLORYS` or `CMEMS` in `config.py`. The supplied configuration uses **GLORYS for initial** and **CMEMS for boundary**.
- Confirm daily coverage, spatial extent, vertical coordinates, velocity orientation, and file names before processing. Script 1 cannot generate *native staggered GLORYS* data; those are separate inputs.

### Tidal atlases: TPXO10 and FES2014

Set `TIDE_MODEL = 'TPXO10'` or `TIDE_MODEL = 'FES2014'`.

**TPXO10:** `TIDE_ATLAS_DIR = DATA_DIR / 'TPXO10'` with the three configured atlas/grid files above. The current constituent list in `config.py` contains 25 names. **Verify that all 25 are actually present in the installed atlas**. For SAL, `TIDE_GENERATE_SAL=True` additionally requires `sal_tpxo10.v2a.nc` at `TIDE_SAL_FILE`. The SAL script may write zero-valued fields for absent constituents; review those results.

**FES2014:** the code expects one NetCDF file per constituent in each of these subdirectories:

```text
Data/FES2014/ocean_tide/m2.nc
Data/FES2014/eastward_velocity/m2.nc
Data/FES2014/northward_velocity/m2.nc
# ... repeat for every configured FES2014 constituent
```

FES2014 is an alternative **boundary-tide** source in script 6. The configured SAL implementation is TPXO10-specific; **do not assume changing `TIDE_MODEL` to FES2014 also changes the SAL atlas**. Disable SAL or supply/verify the intended TPXO10 SAL input and NEMO configuration independently.

Tidal atlases may have access/licensing requirements; obtain them from their official distributors. Inspect amplitude/transport units, harmonic phase conventions, masks and tidal constituent names before running NEMO.

### GloFAS discharge

Script 3 requests monthly files using the configured `GLOFAS_FILENAME`, e.g. `Glofass_global_01_2019.nc`, and checks for `avg_dis` and `valid_time`. Script 9 reads those monthly files, identifies near-coastal outlets, and generates monthly NEMO runoff files with daily records. **The gridded discharge field is not a verified river-mouth inventory**; independently inspect mapped outlet locations and discharge conservation.

### ERA5 and chlorophyll

- Script 2 downloads ERA5 hourly data and produces one daily file **per variable**, including `u10`, `v10`, `t2m`, `d2m`, `msl`, `ssrd`, `strd`, `tp`. Processed files follow the code's `ERA5_1H_<variable>_yYYYYmMMdDD.nc` naming; raw files follow the configurable template with a `raw_` prefix.
- Script 10 requires `CHL_FILE` with finite, separable 2D `nav_lon`/`nav_lat` and a variable named by `CHL_VARIABLE` (default `CHLA`). The code rejects a truly curvilinear chlorophyll grid; do not assume every ocean-colour product is compatible.
- Scripts 8 and 10 prepare SCRIP remapping grids and call compiled NEMO WEIGHTS programs; the Python environment alone is insufficient.

## 4. Configuration: edit `config.py` first

The main settings in the supplied configuration are:

```python
START_DATE = '2019-01-01'
END_DATE = '2019-12-31'
rn_DT = 150.0
NPROCS = 4
INITIAL_TILE_ROWS = 8
INITIAL_SOURCE_TYPE = 'GLORYS'
BOUNDARY_SOURCE_TYPE = 'CMEMS'
TIDE_MODEL = 'TPXO10'
TIDE_GENERATE_BOUNDARY = True
TIDE_GENERATE_SAL = True
TIDE_SAL_TILE_ROWS = 128
CHL_ENABLED = True
NEMO_WEIGHTS_DIR = Path('/path/to/NEMO/tools/WEIGHTS')  # CHANGE THIS
```

**Critical period mismatch in the supplied configuration:** simulation dates cover **2019-01-01 through 2019-12-31**, but CMEMS and ERA5 downloads end **2019-03-31**, while GloFAS downloads cover **2020-10-01 through 2025-12-31**. Change these date ranges or provide existing complete datasets before attempting a full-year run. Script 11 warns about incomplete coverage; that warning does not generate missing data.

`NPROCS` controls script 5 worker processes; `INITIAL_TILE_ROWS` controls row count per initial-condition tile; `TIDE_SAL_TILE_ROWS` controls the full-grid SAL output block size and is not the processor count. With multiprocessing `spawn`, **keep the `if __name__ == '__main__':` entry-point guard in script 5**. Do not remove it: otherwise workers re-import and re-execute the main program, causing recursive worker startup and NetCDF file contention.

## 5. Running all 11 programs

Run each command from the directory containing `config.py`. Steps 1–3 download data and may be run independently/in parallel after setting compatible dates and credentials. The dependency-aware order is shown below.

| Step | Command | Input prerequisites | Main output / purpose |
|---|---|---|---|
| 1 | `python 1_download_cmems.py` | Copernicus Marine credentials, dataset IDs, bounds and dates | Daily CMEMS `gridT/gridS/gridU/gridV/grid2D` NetCDF files under `Data/CMEMS/` |
| 2 | `python 2_download_era5.py` | CDS credentials, ERA5 terms, dates/bounds | Hourly ERA5 data split into daily files under `run/era5/` |
| 3 | `python 3_download_glofass.py` | EWDS credentials, GloFAS access, dates | Monthly `Glofass_global_MM_YYYY.nc` files in `run/glofass_data/` |
| 4 | `python 4_make_boundary_coordinates.py` | `mesh_mask.nc`, configured `RIM_WIDTH` | `run/coordinates.bdy.nc` |
| 5 | `python 5_make_initial.py` | `coordinates.nc`, `mesh_mask.nc`, `domain_cfg.nc`, five first-day GLORYS/CMEMS files | `run/initial_state.nc` |
| 6 | `python 6_make_tides.py` | `coordinates.bdy.nc`, `coordinates.nc`, TPXO/FES atlas, optional TPXO SAL file | Per-constituent `run/tides/*_grd_Z/U/V.nc`; optional full-grid SAL NetCDF |
| 7 | `python 7_make_ocean_boundary.py` | `coordinates.bdy.nc`, grid/geometry files, daily ocean source files covering run | `run/bdy/bdyT_ocean.nc`, `bdyU_ocean.nc`, `bdyV_ocean.nc` |
| 8 | `python 8_make_era5_weights.py` | ERA5 sample, `coordinates.nc`, compiled NEMO WEIGHTS executables | `run/weights/weights_era5_bilin.nc` |
| 9 | `python 9_make_river_runoff.py` | Monthly GloFAS source, `domain_cfg.nc`, `mesh_mask.nc` | `run/runoff/indo36_river_runoff_yYYYYmMM.nc` with daily records |
| 10 | `python 10_make_chlorophyll.py` | Chlorophyll NetCDF, `coordinates.nc`, compiled NEMO WEIGHTS executables | `run/weights/weights_chl_bilin.nc` |
| 11 | `python 11_make_namelist_cfg.py` | `config.py`, **`namelist_cfg_AMM12_template`**, generated input names; monthly runoff files if enabled | `run/namelist_cfg` |

A typical sequential workflow after downloading or staging data:

```bash
conda activate indo36
cd /mnt/d/indo36/Program
python 4_make_boundary_coordinates.py
python 5_make_initial.py
python 6_make_tides.py
python 7_make_ocean_boundary.py
python 8_make_era5_weights.py
python 9_make_river_runoff.py
python 10_make_chlorophyll.py
python 11_make_namelist_cfg.py
```

The scripts write logs under `RUN_DIR/log/`; inspect the corresponding log after each step. The `run/` directory is the preprocessing output location, **not** automatically a complete executable NEMO run directory.

## 6. NEMO namelist template and model setup

Script 11 reads `namelist_cfg_AMM12_template` and modifies its sections, including boundary, tidal, initial, surface forcing, river and chlorophyll settings. **The template is not included in the 12 uploaded Python files.** Obtain the appropriate template from your verified NEMO configuration, place it alongside the scripts, and check that all expected namelist sections exist. The output `run/namelist_cfg` must be reviewed against the **actual NEMO 4.2.2 source and compile options**.

Important model-side items that are **not installed/generated by this environment** include the NEMO binary, required XIOS runtime and XML definitions if your build uses XIOS, executable launch scripts, MPI libraries, a consistent `domain_cfg.nc`, and any other namelists or external forcing demanded by your compiled configuration.

### Outstanding compatibility checks (do not skip)

1. **ERA5 daily file naming and fldread:** the downloader's `ERA5_1H_<var>_yYYYYmMMdDD.nc` files and script 11's generated field descriptors/stems must agree with NEMO's actual time-varying file lookup. Inspect `namelist_cfg` and test a date transition, not just day one.
2. **ERA5 flux units:** script 2 contains conversions for ERA5 accumulated `ssrd`, `strd`, `tp`; inspect resulting units, signs, time semantics and NEMO surface-forcing variable mappings. Confirm the model's snow forcing descriptor refers to a **real variable/file**, not an absent `snow` variable in the precipitation file.
3. **Initial condition:** verify `initial_state.nc` variables, dimensions, masks, vertical coordinates, U/V orientation, and SSH handling against NEMO 4.2.2 initialization conventions.
4. **Ocean boundary:** verify the three boundary files' time coordinates, vector orientation, boundary index/rim conventions, and interpolation at wet points.
5. **Tides/SAL:** verify `z1/z2`, `u1/u2`, `v1/v2` harmonic conventions, constituent selection, SAL variable conventions, and tide-related NEMO switches. In script 6 the SAL output is TPXO10-specific.
6. **River runoff:** ensure source GloFAS dates cover the run and mapped outlets are geographically plausible. The sample `RIVER_SALINITY = 25.0` g/kg is **not freshwater**; change to an appropriate physically justified value if the intended source is river freshwater.
7. **Chlorophyll:** verify source-grid regularity, monthly time axis/climatology metadata, and field-read naming conventions.
8. **Performance:** test script 5 with a small `NPROCS` first. Large `NPROCS` and data stored under WSL `/mnt/d` can cause I/O contention and memory pressure.

These are required validation items; this manual does **not** claim the Python-generated files have already passed a NEMO integration test.

## 7. Input checks and troubleshooting

Inspect source files before execution:

```bash
ncdump -h /mnt/d/indo36/Data/coordinates.nc | head -80
ncdump -h /mnt/d/indo36/Data/mesh_mask.nc | head -80
ncdump -h /mnt/d/indo36/Data/domain_cfg.nc | head -80
ncdump -h /mnt/d/indo36/Data/TPXO10/h_tpxo10.v2a.nc | head -80
ncdump -h /mnt/d/indo36/Data/TPXO10/u_tpxo10.v2a.nc | head -80
ncdump -h /mnt/d/indo36/Data/TPXO10/grid_tpxo10v2a.nc | head -80
```

Quick Python checks:

```bash
python - <<'PY'
from pathlib import Path
import config as C
for name in ('COORDINATES', 'MESH_MASK', 'DOMAIN_CFG'):
    p = Path(getattr(C, name))
    print(f'{name}: {p} exists={p.is_file()}')
print('Initial source:', C.INITIAL_SOURCE_TYPE)
print('Boundary source:', C.BOUNDARY_SOURCE_TYPE)
print('Tide model:', C.TIDE_MODEL)
print('Period:', C.START_DATE, C.END_DATE)
PY
```

**`PermissionError` / repeated `START workers=` during script 5:** check that its multiprocessing entry-point guard is present. Never remove it when using `ProcessPoolExecutor(..., mp_context=get_context('spawn'))`. Stop failed worker processes and delete only a confirmed stale `.pending.nc` file before retrying.

**`FileNotFoundError`:** compare the exact logged path with `config.py`, including filename date conventions and the correct ocean source type.

**`Missing NEMO WEIGHTS executable`:** compile `scripgrid.exe`, `scrip.exe`, `scripshape.exe` and set `NEMO_WEIGHTS_DIR` to their actual directory.

**`Tide target outside source rectangle`:** inspect atlas longitude convention, geographic coverage, and selected domain. Do not suppress this error without verifying the interpolation geometry.

**`No ERA5 daily files found`:** script 8 searches `ERA5_PREPARED_DIR` for either raw daily files or processed files matching `ERA5_1H_<var>_y????m??d??.nc`.

**`Missing template namelist section`:** supply the intended `namelist_cfg_AMM12_template` with the expected section names; do not replace it with an arbitrary empty file.

## 8. Data sources and external documentation

- [NEMO ocean model](https://www.nemo-ocean.eu/) — NEMO 4.2.2 source/build and namelist documentation.
- [Copernicus Marine](https://marine.copernicus.eu/) — CMEMS daily physical ocean products and credentials.
- [Climate Data Store](https://cds.climate.copernicus.eu/) — ERA5 hourly reanalysis.
- [Copernicus Emergency Management / EWDS](https://ewds.climate.copernicus.eu/) — GloFAS historical discharge products.
- [TPXO tidal models](https://www.tpxo.net/) — tidal atlas distribution, constituent availability, units and access terms.
- [FES tides](https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html) — FES atlas information and distribution.
- [ESA Ocean Colour CCI](https://www.oceancolour.org/) — chlorophyll source products.

## 9. Publication and reproducibility

Record the exact source dataset versions, download dates, NEMO revision, `config.py`, compiler and MPI build options, tide atlas version, and checksums of the produced NetCDF files. Keep large source and output NetCDF files outside the Git repository. Include licensing/attribution required by the original data providers when publishing outputs.

---

**Author:** Dr Dwiyoga Nugroho — Toulouse, November 2026.
