# Cross-Matching Gaia

## ADQL at Gaia Archive

### Mathur17

1. Generate VOTables and upload to Gaia Archive
2. Log into Gaia archive as `epetigur`
3. Perform the following joins that query all gaia sources within 8 arcsec of each CKS source, or stellar17 source.
4. Give the job a sensible name like `xmatch_cks_gaiadr2`
5. Run the query
6. Download results as a votable
7. Move into data

Same as CKS, but with this query

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

Run with top 1000 to verify that it works 

Should take about 3min to run 