# Cross-Matching Gaia

## ADQL at Gaia Archive

### Mathur17

1. Save following keys to a csv file
   id_kic
   id_tmass
   ra
   dec

2. Log into Gaia archive as `epetigur`



3. Perform the cross match following joins that query all gaia sources within 8
   arcsec of each CKS source, or stellar17 source.

Merge based on distance

```
SELECT *,distance(
  POINT('ICRS', m17.m17_ra, m17.m17_dec),
  POINT('ICRS', gaia.ra, gaia.dec)) AS dist
FROM gaiadr1.gaia_source AS gaia, 
     user_epetigur.m17 AS m17, 
     gaiadr2.ruwe as ruwe
WHERE 1=CONTAINS(
  POINT('ICRS', m17.m17_ra, m17.m17_dec),
  CIRCLE('ICRS', gaia.ra, gaia.dec, 0.00222)
)
```

Merge based on gaia-tmass-kic crossmatch

Run with top 1000 to verify that it works 

Should take about 3min to run 

```
SELECT TOP 10
m17.id_kic as id_kic,
m17.id_tmass as id_tmass,
gaia.source_id as id_gaia2,
gaia.ra as gaia2_ra, 
gaia.dec as gaia2_dec,
gaia.parallax as gaia2_sparallax, 
gaia.parallax_error as gaia2_sparallax_err, 
gaia.phot_g_mean_flux as gaia2_gflux,
gaia.phot_g_mean_flux_error as gflux_err,
gaia.phot_g_mean_mag as gaia2_gmag,
gaia.phot_bp_mean_flux as gaia2_bpflux,
gaia.phot_bp_mean_flux_error as gaia2_bpflux_err,
gaia.phot_bp_mean_mag as gaia2_bpmag,
gaia.phot_rp_mean_flux as gaia2_rpflux,
gaia.phot_rp_mean_flux_error as gaia2_rpflux_err,
gaia.phot_rp_mean_mag as gaia2_rpmag,
gaia.parallax_over_error as gaia2_sparallax_over_err,
gaia.astrometric_excess_noise as gaia2_astrometric_excess_noise,
gaia.astrometric_excess_noise_sig as gaia2_astrometric_excess_noise_sig,
gaia.teff_val as gaia2_steff,
gaia.teff_percentile_upper - gaia.teff_val as gaia2_steff_err1,
gaia.teff_percentile_lower - gaia.teff_val as gaia2_steff_err2,
gaia.radius_val as gaia2_steff,
gaia.radius_percentile_upper - gaia.radius_val as gaia2_steff_err1,
gaia.radius_percentile_lower - gaia.radius_val as gaia2_steff_err2,
ruwe.ruwe as gaia2_ruwe,
xmatch.number_of_neighbours as xm_number_of_neighbours,
xmatch.number_of_mates as xm_number_of_mates,
xmatch.best_neighbour_multiplicity as xm_best_neighbour_multiplicity
xmatch.best_neighbour_multiplicity as xm_best_neighbour_multiplicity
FROM gaiadr2.gaia_source AS gaia
INNER JOIN gaiadr2.ruwe AS ruwe 
ON gaia.source_id = ruwe.source_id
INNER JOIN gaiadr2.tmass_best_neighbour AS xmatch
ON gaia.source_id = xmatch.source_id
INNER JOIN gaiadr1.tmass_original_valid AS tmass
ON tmass.tmass_oid = xmatch.tmass_oid
INNER JOIN user_epetigur.m17_tmass as m17
ON tmass.designation=m17.id_tmass
INNER JOIN external.gaiadr2_geometric_distance as bj
ON gaia.source_id = bj.source_id
```

4. Give the job a sensible name like `xmatch_cks_gaiadr2`
5. Run the query
6. Download results as votable (because they zip it server side)
7. Move into data




