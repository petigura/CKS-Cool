# CKS-Cool

## Use conda environment

source activate ckscool

## Dependencies

### had to install the following packages with conda

conda install scipy
conda install matplotlib
conda install pandas
conda install astropy
conda install h5py
pip install lmfit
pip install ebfpy
conda install -c conda-forge healpy # needed for dustmaps also got healpy==1.11 to work
pip install dustmaps 
pip install mwdust # needed for dustmaps
pip install pyephem
conda install seaborn
conda install scikit-learn
conda install pytables

# Copy over the table cache from CKS-Gaia

$ cp ../CKS-Gaia/load_table_cache.hdf data/cksgaia_cache.hdf

Repo for the CKS-Cool project includes code to

- [Construct target list](docs/observing.md)
- [Plan/execute observations](docs/observing.md)
- Compute stellar/planet parameters
- Access the CKS-Cool dataset

## Compute Stellar/Planetary Properties

1. Download smemp and smsyn catalogs to `data/`

```
rsync -av cadence:/data/user/petigura/public_html/smemp/specmatch_results.csv data/specmatch-emp_results.csv
rsync -av cadence:/data/user/petigura/public_html/smsyn/specmatch_results.csv data/specmatch-syn_results.csv
```

2. Run the isoclassify code

First create the batch processing files. In total there are about 900 stars that pass the photometric only cuts.

```
run_ckscool.py create-iso-batch 
source bin/create_tot.sh
```

Then run them in parallel. Running isoclassify takes about 3s per star / core. Can process in about 7*3 min with six cores.

### Test first 9 from each method

```
head `ls isoclassify*tot` | grep mkdir | parallel
```

### Run them in batch with

cat isoclassify*tot | grep mkdir | parallel



```
run_ckscool.py create-iso-table
```

3. Generate ReaMatch table

Run `ReaMatch.ipynb` and copy output file to `reamatch.csv` to `~/Dropbox/CKS-Cool/hires/reamatch.csv`. Howard Isaacson will then add the appropriate RM designations to file.

4. Generate representative spectra figure

Run the 3_Spectra-Figure ipython notebook

```
rsync -av --progress --files-from=data/fig_spectra/fig_spectra-files.txt cadence:/ data/fig_spectra/ 
```

## Access CKS-Cool dataset

The full list of star and planet properties are in 

`data/ckscool-planets-cuts.csv `

see

`data/column-definitions.txt`

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




