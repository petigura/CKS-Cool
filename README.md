# CKS-Cool

Repo for planning observations and initial catalog selection

# Generate HIRES scripts

## Update progress during run. 

cps sync # update logs
run_ckscool.py create-val stat

- Download all the NIRC2 ascii files from KOA with Petigura as PI. 
- save them in data/nirc2/ascii


# 

# Catalog paper

1. Download smemp and smsyn catalogs to data/


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

run_ckscool.py create-iso-table

3. Generate ReaMatch table

Run ReaMatch ipython notebook. Copy ~/Dropbox/CKS-Cool/hires/reamatch.csv

4. Download representative spectra

Run the 3_Spectra-Figure ipython notebook

rsync -av --progress --files-from=data/fig_spectra/fig_spectra-files.txt cadence:/ data/fig_spectra/ 