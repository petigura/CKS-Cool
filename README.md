# CKS-Cool

## Computing environment

Use the ckscool environment on Petigura's laptop

```bash
conda activate ckscool 
```
conda install numpy==1.15.4 # this avoids the ValueError: cannot set WRITEABLE flag to True of this array #24839
conda install scipy matplotlib astropy pandas d seaborn scikit-learn pytables
conda install -c conda-forge healpy # needed for dustmaps also got healpy==1.11 to work
conda install joblib # needed for occurrence valley work
pip install h5py # conda doesn't seem to work
pip install mwdust 
pip install pyephem
pip install lmfit
pip install ebfpy

isoclassify on f6f16e

## Cookbook for target list construction 

- [Construct target list](docs/observing.md)
- [Plan/execute observations](docs/observing.md)

## Cookbook for occurrence analysis

- Compute stellar/planet parameters
- Access the CKS-Cool dataset

# Copy over the table cache from CKS-Gaia

$ cp ../CKS-Gaia/load_table_cache.hdf data/cksgaia_cache.hdf


## Compute Stellar/Planetary Properties

1. Download smemp and smsyn catalogs to `data/`

https://jump.caltech.edu/explorer/58/

2. Run the isoclassify code on combined DR1+DR2 dataset.

Run the isoclassify code in three modes on three different
spectroscopic parameters. On the combined DR1+DR2 data

First generate the csv files of stellar parameters

```
bin/run_ckscool.py create-iso-batch # creates 9 csv files
```

### Test first 9 from each method

```bash
isoclassify multiproc direct -1 data/isoclassify-cks1-direct.csv.16 isoclassify/cks1/direct.csv --baseoutdir isoclassify/cks1/direct/
isoclassify multiproc grid -1 data/isoclassify-cks1-grid-parallax-yes.csv.16 isoclassify/cks1/grid-parallax-yes.csv --baseoutdir isoclassify/cks1/grid-parallax-yes/
isoclassify multiproc grid -1 data/isoclassify-cks1-grid-parallax-no.csv.16 isoclassify/cks1/grid-parallax-no.csv --baseoutdir isoclassify/cks1/grid-parallax-no/
isoclassify multiproc direct -1 data/isoclassify-smsyn-direct.csv.16 isoclassify/smsyn/direct.csv --baseoutdir isoclassify/smsyn/direct/
isoclassify multiproc grid -1 data/isoclassify-smsyn-grid-parallax-yes.csv.16 isoclassify/smsyn/grid-parallax-yes.csv --baseoutdir isoclassify/smsyn/grid-parallax-yes/
isoclassify multiproc grid -1 data/isoclassify-smsyn-grid-parallax-no.csv.16 isoclassify/smsyn/grid-parallax-no.csv --baseoutdir isoclassify/smsyn/grid-parallax-no/
isoclassify multiproc direct -1 data/isoclassify-smemp-direct.csv.16 isoclassify/smemp/direct.csv --baseoutdir isoclassify/smemp/direct/
isoclassify multiproc grid -1 data/isoclassify-smemp-grid-parallax-yes.csv.16 isoclassify/smemp/grid-parallax-yes.csv --baseoutdir isoclassify/smemp/grid-parallax-yes/
isoclassify multiproc grid -1 data/isoclassify-smemp-grid-parallax-no.csv.16 isoclassify/smemp/grid-parallax-no.csv --baseoutdir isoclassify/smemp/grid-parallax-no/
```

Look at the logs and confirm isoclassify is behaving right. Run them on Erik's laptop.


```bash
isoclassify multiproc direct 6 data/isoclassify-cks1-direct.csv isoclassify/cks1/direct.csv --baseoutdir isoclassify/cks1/direct/ --plot none
isoclassify multiproc grid 6 data/isoclassify-cks1-grid-parallax-yes.csv isoclassify/cks1/grid-parallax-yes.csv --baseoutdir isoclassify/cks1/grid-parallax-yes/ --plot none
isoclassify multiproc grid 6 data/isoclassify-cks1-grid-parallax-no.csv isoclassify/cks1/grid-parallax-no.csv --baseoutdir isoclassify/cks1/grid-parallax-no/ --plot none
isoclassify multiproc direct 6 data/isoclassify-smsyn-direct.csv isoclassify/smsyn/direct.csv --baseoutdir isoclassify/smsyn/direct/  --plot none
isoclassify multiproc grid 6 data/isoclassify-smsyn-grid-parallax-yes.csv isoclassify/smsyn/grid-parallax-yes.csv --baseoutdir isoclassify/smsyn/grid-parallax-yes/ --plot none
isoclassify multiproc grid 6 data/isoclassify-smsyn-grid-parallax-no.csv isoclassify/smsyn/grid-parallax-no.csv --baseoutdir isoclassify/smsyn/grid-parallax-no/ --plot none
isoclassify multiproc direct 6 data/isoclassify-smemp-direct.csv isoclassify/smemp/direct.csv --baseoutdir isoclassify/smemp/direct/ --plot none
isoclassify multiproc grid 6 data/isoclassify-smemp-grid-parallax-yes.csv isoclassify/smemp/grid-parallax-yes.csv --baseoutdir isoclassify/smemp/grid-parallax-yes/ --plot none
isoclassify multiproc grid 6 data/isoclassify-smemp-grid-parallax-no.csv isoclassify/smemp/grid-parallax-no.csv --baseoutdir isoclassify/smemp/grid-parallax-no/ --plot none
```

The create isoclassify tables

```bash
run_ckscool.py create-iso-table
```


3. Generate ReaMatch table

Run `ReaMatch.ipynb` and copy output file to `reamatch.csv` to `~/Dropbox/CKS-Cool/hires/reamatch.csv`. Howard Isaacson will then add the appropriate RM designations to file.


4. Create final table of stellar / planet parameters





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






# Old notes

## Cross-match the CKS stars with Gaia IDs. This is now accomplished using the Berger designatinos

Follow instructions [here](docs/gaia-xmatch.md)



# Notes for running isoclassify on cadence

Note when running on cadence, there was a really weird issue with h5py. Where it was taking 30s to read in the Combined Dustmap
