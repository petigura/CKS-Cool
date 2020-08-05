# CKS-Cool

## Use conda environment

source activate ckscool

## Dependencies

### had to install the following packages with conda

conda install numpy==1.15.4 # this avoids the ValueError: cannot set WRITEABLE flag to True of this array #24839
conda install scipy matplotlib astropy pandas d seaborn scikit-learn pytables
conda install -c conda-forge healpy # needed for dustmaps also got healpy==1.11 to work
conda install joblib
#
pip install h5py # conda doesn't seem to work
pip install mwdust 
pip install pyephem
pip install lmfit
pip install ebfpy

isoclassify on f6f16ef1f90c57893268f6f9d6da4fddb00142ec

Note when running on cadence, there was a really weird issue with h5py. Where it was taking 30s to read in the Combined Dustmap

# Copy over the table cache from CKS-Gaia

$ cp ../CKS-Gaia/load_table_cache.hdf data/cksgaia_cache.hdf

Repo for the CKS-Cool project includes code to

- [Construct target list](docs/observing.md)
- [Plan/execute observations](docs/observing.md)
- Compute stellar/planet parameters
- Access the CKS-Cool dataset

## Cross-match the CKS stars with Gaia IDs.

Follow instructions [here](docs/gaia-xmatch.md)

## Compute Stellar/Planetary Properties

1. Download smemp and smsyn catalogs to `data/`

```
rsync -av cadence:/data/user/doppler/public_html/smemp/specmatch_results.csv data/specmatch-emp_results.csv
rsync -av cadence:/data/user/doppler/public_html/smsyn/specmatch_results.csv data/specmatch-syn_results.csv
```


2. Create kbc table

A list of all Kepler templates from HIRES

```
run_ckscool.py create-kbc
```

2. Run the isoclassify code

First create the batch processing files. In total there are about 900 stars that pass the photometric only cuts.

export DUST_DIR=/data/user/petigura/dustdir/
```
run_ckscool.py create-iso-batch 
```

```
python bin/create_tot.py
```

Then run them in parallel. Running isoclassify takes about 3s per star / core. Can process in about 7*3 min with six cores.

### Test first 9 from each method

```
head `ls isoclassify*tot` | grep mkdir | parallel
```

### Run them in batch on cadence 

Create the tot files. See note about h5py. 

```
source bin/create_tot.sh
```

Run isoclassify on a screen session. If I run with all the cores on
cadence I get a seg fault.

```
cat isoclassify*tot | grep mkdir | parallel -j 48 
```

I also need to set this environment variable or else I get a resource unavailable error from hdf

export HDF5_USE_FILE_LOCKING=FALSE 
DUST_DIR=/data/user/petigura/dustdir/

Notes:

I tried to use the bayestar interface, but I got a I get a "too many
requests" error 


### Monitor job progress with

```
 echo "number of log files"; find isoclassify/ -name "*.log" | wc -l  ; echo "number of csv files"; grep csv `find isoclassify/ -name "*.log" ` | grep created | wc -l 
```

First number is number of log files, second numbers is how many csv files created.

```
run_ckscool.py create-iso-table
```

Run it on cadence

3. Generate ReaMatch table

Run `ReaMatch.ipynb` and copy output file to `reamatch.csv` to `~/Dropbox/CKS-Cool/hires/reamatch.csv`. Howard Isaacson will then add the appropriate RM designations to file.

4. Generate representative spectra figure

Run the 3_Spectra-Figure ipython notebook

```
rsync -av --progress --files-from=data/fig_spectra/fig_spectra-files.txt cadence:/ data/fig_spectra/ 
```

## Access CKS-Cool dataset

The full list of star and planet properties are in `data/ckscool-planets-cuts.csv ` see `data/column-definitions.txt` for a description of the columns.

The `is*` columns correspond to cuts. See the ckscool/cuts.py for additional info. The radius gap be comes more clear when one adopts a koi_impact of > 0.7 or 0.8.




For a description of the columns

## Generate HDF version of DR25 chains.

First run

```
run_ckscool.py create-chain-hdf
```

To read the chain info from its ascii format

Notes
- Chains take about 2s to read and store in hdf format or 4 hours for all 8000
- 4MB per chain or 32 GB for all 8000.
- There are roughly ~100 KOIs for which there are no chains

```
run_ckscool.py create-chain-summary # stores the precentile summary
```

Notes 
- Takes about 30 min to complete.
