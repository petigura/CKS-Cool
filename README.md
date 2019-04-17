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
conda install -c conda-forge healpy # needed for dustmaps
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

```
run_ckscool.py create-iso-batch 
isoclassify batch direct data/isoclassify-direct.csv -o isoclassify/direct/ > isoclassify-direct.tot
isoclassify batch grid data/isoclassify-grid-parallax-yes.csv  -o isoclassify/grid-parallax-yes/ > isoclassify-grid-parallax-yes.tot
isoclassify batch grid data/isoclassify-grid-parallax-no.csv -o isoclassify/grid-parallax-no/ > isoclassify-grid-parallax-no.tot

cat isoclassify-direct.tot | parallel -j 6 
cat isoclassify-grid-parallax-no.tot | parallel -j 6  # takes about 20min on Erik's laptop
cat isoclassify-grid-parallax-yes.tot | parallel -j 6 # takes about 20min on Erik's laptop
```

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




