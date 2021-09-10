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


### Assemble the data

1. Copy over CKS-Gaia HDF file

```bash
cp ../CKS-Gaia/load_table_cache.hdf data/cksgaia_cache.hdf
```

2. Generate a symlink to `dr25-chains_trimmed-thinned.hdf`

Note the DR25 chains are created in `Kepler-Radius-Ratio` repo. See
its readme file. Users looking to rerun code, will need to get the
datafile from Petigura.

### Compute Stellar/Planetary Properties

1. Download smemp and smsyn catalogs to `data/CKS_Spectroscopic_Parameters.csv`

```
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

conda deactivate

The create isoclassify tables

```bash
run_ckscool.py create-iso-table
```


3. Generate ReaMatch table

Run `ReaMatch.ipynb` and copy output file to `reamatch.csv` to `~/Dropbox/CKS-Cool/hires/reamatch.csv`. Howard Isaacson will then add the appropriate RM designations to file.


4. Create final table of stellar / planet parameters

5. Generate representative spectra figure

Run the 3_Spectra-Figure ipython notebook

```
rsync -av --progress --files-from=data/fig_spectra/fig_spectra-files.txt cadence:/ data/fig_spectra/ 
```



## Create plots and build paper

## Run gapfitting jupyter to generate gap fits. 


## Access CKS-Cool dataset

The full list of star and planet properties are in `data/ckscool-planets-cuts.csv ` see `data/column-definitions.txt` for a description of the columns.

The `is*` columns correspond to cuts. See the ckscool/cuts.py for additional info.

## Other notes

1. Running isoclassify on cadence

Note when running on cadence, there was a really weird issue with
h5py. Where it was taking 30s to read in the Combined Dustmap

## Download all Kepler headers (this is needed for the dilution cut)


List of wget scripts from all quarters

https://exoplanetarchive.ipac.caltech.edu/bulk_data_download/Kepler_Quarterly_wget.tar.gz


## Combine wget scripts and cat them into one file and only pull the fits files.

cat ~/Downloads/Kepler_Quarterly_wget/Kepler_Q* | grep fits > Kepler_wget.bat

run this script on cadence
 
Downloads about 150 lightcurves/s. Should take about 4 hours to download all.

##

scrape_headers.py

CKS-Cool xmatch_gaia2_m17_ruwe_tmass

Here's the Gaia query.

SELECT m17.id_kic as id_kic, m17.id_tmass as id_tmass, gaia.source_id as id_gaia2, gaia.ra as gaia2_ra, gaia.dec as gaia2_dec, gaia.parallax as gaia2_sparallax, gaia.parallax_error as gaia2_sparallax_err, gaia.phot_g_mean_flux as gaia2_gflux, gaia.phot_g_mean_flux_error as gflux_err, gaia.phot_g_mean_mag as gaia2_gmag, gaia.phot_bp_mean_flux as gaia2_bpflux, gaia.phot_bp_mean_flux_error as gaia2_bpflux_err, gaia.phot_bp_mean_mag as gaia2_bpmag, gaia.phot_rp_mean_flux as gaia2_rpflux, gaia.phot_rp_mean_flux_error as gaia2_rpflux_err, gaia.phot_rp_mean_mag as gaia2_rpmag, gaia.parallax_over_error as gaia2_sparallax_over_err, gaia.astrometric_excess_noise as gaia2_astrometric_excess_noise, gaia.astrometric_excess_noise_sig as gaia2_astrometric_excess_noise_sig, gaia.teff_val as gaia2_steff, gaia.teff_percentile_upper - gaia.teff_val as gaia2_steff_err1, gaia.teff_percentile_lower - gaia.teff_val as gaia2_steff_err2, gaia.radius_val as gaia2_srad, gaia.radius_percentile_upper - gaia.radius_val as gaia2_srad_err1, gaia.radius_percentile_lower - gaia.radius_val as gaia2_srad_err2, ruwe.ruwe as gaia2_ruwe, xmatch.number_of_neighbours as xm_number_of_neighbours, xmatch.number_of_mates as xm_number_of_mates, xmatch.best_neighbour_multiplicity as xm_best_neighbour_multiplicity, bj.r_lo as bj_dist, gaia.radius_percentile_upper - gaia.radius_val as gaia2_steff_err1, gaia.radius_percentile_lower - gaia.radius_val as gaia2_steff_err2, bj.r_hi - bj.r_est as bj_dist_err1, bj.r_lo - bj.r_est as bj_dist_err2 FROM gaiadr2.gaia_source AS gaia INNER JOIN gaiadr2.ruwe AS ruwe ON gaia.source_id = ruwe.source_id INNER JOIN gaiadr2.tmass_best_neighbour AS xmatch ON gaia.source_id = xmatch.source_id INNER JOIN gaiadr1.tmass_original_valid AS tmass ON tmass.tmass_oid = xmatch.tmass_oid INNER JOIN user_epetigur.m17_tmass as m17 ON tmass.designation=m17.id_tmass INNER JOIN external.gaiadr2_geometric_distance as bj ON gaia.source_id = bj.source_id

## return the same gaia columns for all sources within 10 arcsec of target

SELECT gaia.source_id as id_gaia2, gaia.ra as gaia2_ra, gaia.dec as gaia2_dec, gaia.parallax as gaia2_sparallax, gaia.parallax_error as gaia2_sparallax_err, gaia.phot_g_mean_flux as gaia2_gflux, gaia.phot_g_mean_flux_error as gflux_err, gaia.phot_g_mean_mag as gaia2_gmag, gaia.phot_bp_mean_flux as gaia2_bpflux, gaia.phot_bp_mean_flux_error as gaia2_bpflux_err, gaia.phot_bp_mean_mag as gaia2_bpmag, gaia.phot_rp_mean_flux as gaia2_rpflux, gaia.phot_rp_mean_flux_error as gaia2_rpflux_err, gaia.phot_rp_mean_mag as gaia2_rpmag, gaia.parallax_over_error as gaia2_sparallax_over_err, gaia.astrometric_excess_noise as gaia2_astrometric_excess_noise, gaia.astrometric_excess_noise_sig as gaia2_astrometric_excess_noise_sig, gaia.teff_val as gaia2_steff, gaia.teff_percentile_upper - gaia.teff_val as gaia2_steff_err1, gaia.teff_percentile_lower - gaia.teff_val as gaia2_steff_err2, gaia.radius_val as gaia2_srad, gaia.radius_percentile_upper - gaia.radius_val as gaia2_srad_err1, gaia.radius_percentile_lower - gaia.radius_val as gaia2_srad_err2,
distance(
  POINT('ICRS', subquery.gaia2_ra, subquery.gaia2_dec),
  POINT('ICRS', gaia.ra, gaia.dec)) AS dist
FROM gaiadr2.gaia_source AS gaia
JOIN
(
SELECT gaia.ra as gaia2_ra, gaia.dec as gaia2_dec
FROM gaiadr2.gaia_source AS gaia INNER JOIN gaiadr2.ruwe AS ruwe ON gaia.source_id = ruwe.source_id INNER JOIN gaiadr2.tmass_best_neighbour AS xmatch ON gaia.source_id = xmatch.source_id INNER JOIN gaiadr1.tmass_original_valid AS tmass ON tmass.tmass_oid = xmatch.tmass_oid INNER JOIN user_epetigur.m17_tmass as m17 ON tmass.designation=m17.id_tmass INNER JOIN external.gaiadr2_geometric_distance as bj ON gaia.source_id = bj.source_id
) AS subquery
ON
1=CONTAINS(
  POINT('ICRS', subquery.gaia2_ra, subquery.gaia2_dec),
  CIRCLE('ICRS', gaia.ra, gaia.dec, 0.002777777778))







SELECT subquery.id_gaia2 as id_gaia2_primary, gaia.source_id as id_gaia2, gaia.ra as gaia2_ra, gaia.dec as gaia2_dec, gaia.parallax as gaia2_sparallax, gaia.parallax_error as gaia2_sparallax_err, gaia.phot_g_mean_flux as gaia2_gflux, gaia.phot_g_mean_flux_error as gflux_err, gaia.phot_g_mean_mag as gaia2_gmag, gaia.phot_bp_mean_flux as gaia2_bpflux, gaia.phot_bp_mean_flux_error as gaia2_bpflux_err, gaia.phot_bp_mean_mag as gaia2_bpmag, gaia.phot_rp_mean_flux as gaia2_rpflux, gaia.phot_rp_mean_flux_error as gaia2_rpflux_err, gaia.phot_rp_mean_mag as gaia2_rpmag, gaia.parallax_over_error as gaia2_sparallax_over_err, gaia.astrometric_excess_noise as gaia2_astrometric_excess_noise, gaia.astrometric_excess_noise_sig as gaia2_astrometric_excess_noise_sig, gaia.teff_val as gaia2_steff, gaia.teff_percentile_upper - gaia.teff_val as gaia2_steff_err1, gaia.teff_percentile_lower - gaia.teff_val as gaia2_steff_err2, gaia.radius_val as gaia2_srad, gaia.radius_percentile_upper - gaia.radius_val as gaia2_srad_err1, gaia.radius_percentile_lower - gaia.radius_val as gaia2_srad_err2,
distance(
  POINT('ICRS', subquery.gaia2_ra, subquery.gaia2_dec),
  POINT('ICRS', gaia.ra, gaia.dec)) AS dist
FROM gaiadr2.gaia_source AS gaia
JOIN
(
SELECT TOP 10 gaia.source_id as id_gaia2, gaia.ra as gaia2_ra, gaia.dec as gaia2_dec
FROM gaiadr2.gaia_source AS gaia INNER JOIN gaiadr2.ruwe AS ruwe ON gaia.source_id = ruwe.source_id INNER JOIN gaiadr2.tmass_best_neighbour AS xmatch ON gaia.source_id = xmatch.source_id INNER JOIN gaiadr1.tmass_original_valid AS tmass ON tmass.tmass_oid = xmatch.tmass_oid INNER JOIN user_epetigur.m17_tmass as m17 ON tmass.designation=m17.id_tmass INNER JOIN external.gaiadr2_geometric_distance as bj ON gaia.source_id = bj.source_id
) AS subquery
ON
1=CONTAINS(
  POINT('ICRS', subquery.gaia2_ra, subquery.gaia2_dec),
  CIRCLE('ICRS', gaia.ra, gaia.dec, 0.002777777778))
