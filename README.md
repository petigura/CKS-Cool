# CKS-Cool

1. [Computing Environemnt](#computing)
2. [Target List Construction](#targetlist)
3. [Cookbook For CKS-X](#cksx)
4. [Data Access](#dataaccess)
5. [Isoclassify Cadence](#isoclassify-cadence)
6. [Kepler Headers](#kepler-headers)


## Computing environment <a name="computing"></a>

Use the ckscool environment on Petigura's laptop

```bash
conda activate ckscool 
conda install numpy==1.15.4 # this avoids the ValueError: cannot set WRITEABLE flag to True of this array #24839
conda install scipy matplotlib astropy pandas d seaborn scikit-learn pytables
conda install -c conda-forge healpy # needed for dustmaps also got healpy==1.11 to work
conda install joblib # needed for occurrence valley work
pip install h5py # conda doesn't seem to work
pip install mwdust 
pip install pyephem
pip install lmfit
pip install ebfpy
```

isoclassify on f6f16e

## Cookbook for target list construction <a name="targetlist"></a>

1. [Construct target list](docs/observing.md)
1. [Plan/execute observations](docs/observing.md)

## Cookbook for occurrence analysis <a name="cksx"></a>

- Compute stellar/planet parameters
- Access the CKS-Cool dataset


### Assemble the data

1. Copy over CKS-Gaia HDF file

   ```bash
   cp ../CKS-Gaia/load_table_cache.hdf data/cksgaia_cache.hdf
   ```

2. Generate a symlink to `dr25-chains_trimmed-thinned.hdf`

   Note the DR25 chains are created in `Kepler-Radius-Ratio` repo. See
   its readme file. Users looking to rerun code, will need to get the
   datafile from Petigura.

3. Query Gaia parameters. Instructions [here](docs/gaia-xmatch.md)

### Compute Stellar/Planetary Properties

1. Download smemp and smsyn catalogs to `data/CKS_Spectroscopic_Parameters.csv`

   ```bash
   https://jump.caltech.edu/explorer/58/
   mv ~/Downloads/CKS_Spectroscopic_Parameters.csv data/
   ```

2. Run the isoclassify code on combined DR1+DR2 dataset.

   Run the isoclassify code in three modes on three different
   spectroscopic parameters. On the combined DR1+DR2 data

   First generate the csv files of stellar parameters

   ```
   bin/run_ckscool.py create-iso-batch # creates 9 csv files
   ```

   Note, if you need to run just a few new stars, add line in`run_ckscool.py`

   conda activate isoclassify

   Look at the logs and confirm isoclassify is behaving right. Run them on Erik's laptop.

   ```bash
   isoclassify multiproc direct 6 data/isoclassify-smsyn-direct.csv isoclassify/smsyn/direct.csv --baseoutdir isoclassify/smsyn/direct/  --plot none
   isoclassify multiproc grid 6 data/isoclassify-smsyn-grid-parallax-yes.csv isoclassify/smsyn/grid-parallax-yes.csv --baseoutdir isoclassify/smsyn/grid-parallax-yes/ --plot none
   isoclassify multiproc grid 6 data/isoclassify-smsyn-grid-parallax-no.csv isoclassify/smsyn/grid-parallax-no.csv --baseoutdir isoclassify/smsyn/grid-parallax-no/ --plot none
   isoclassify multiproc direct 6 data/isoclassify-smemp-direct.csv isoclassify/smemp/direct.csv --baseoutdir isoclassify/smemp/direct/ --plot none
   isoclassify multiproc grid 6 data/isoclassify-smemp-grid-parallax-yes.csv isoclassify/smemp/grid-parallax-yes.csv --baseoutdir isoclassify/smemp/grid-parallax-yes/ --plot none
   isoclassify multiproc grid 6 data/isoclassify-smemp-grid-parallax-no.csv isoclassify/smemp/grid-parallax-no.csv --baseoutdir isoclassify/smemp/grid-parallax-no/ --plot none
   ```

   ```bash
   conda deactivate
   ```

   Them create isoclassify tables

   ```bash
   run_ckscool.py create-iso-table
   ```


### Generate ReaMatch table

Run `ReaMatch.ipynb` and copy output file to `reamatch.csv` to `~/Dropbox/CKS-Cool/hires/reamatch.csv`. Howard Isaacson will then add the appropriate RM designations to file.


### Generate representative spectra figure

   Run the 3_Spectra-Figure ipython notebook

   ```
   rsync -av --progress --files-from=data/fig_spectra/fig_spectra-files.txt cadence:/ data/fig_spectra/ 
   ```
   
### Run the gap fitting code.

   Takes about 1 hour to run with 1000 samples

   ```bash 
   bin/run_ckscool.py create-gradient # stores results in analysis/ to do a fresh run, remove or move this directory
   ```

   Calling this will generate a bunch of other cached results

### Create plots and build paper

run_ckscool.py build plot all 
run_ckscool.py build table all 
run_ckscool.py build val all 
run_ckscool.py build csv all 

http://localhost:8889/notebooks/gapfitting-comparison/gapfitting.ipynb


## Access to the CKS-Cool dataset <a name="dataaccess"></a>

The full list of star and planet properties are in `data/ckscool-planets-cuts.csv ` see `data/column-definitions.txt` for a description of the columns.

The `is*` columns correspond to cuts. See the ckscool/cuts.py for additional info.


## Running isoclassify on cadence <a name="isoclassify-cadence"></a>

Note when running on cadence, there was a really weird issue with
h5py. Where it was taking 30s to read in the Combined Dustmap

## Download all Kepler headers (this is needed for the dilution cut) <a name="kepler-headers"></a>

1. Get bactch download script [here](https://exoplanetarchive.ipac.caltech.edu/bulk_data_download/Kepler_Quarterly_wget.tar.gz)
2. Combine wget scripts and cat them into one file and only pull the fits files.

   ```bash
   cat ~/Downloads/Kepler_Quarterly_wget/Kepler_Q* | grep fits > Kepler_wget.bat
   ```
3. Download onto `cadence:/data/user/petigura/lightcurve-all` cadence Should take about 4 hours to download all.

   ```
   scrape_headers.py
   ```
