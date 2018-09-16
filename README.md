# CKS-Cool

Repo for planning observations and initial catalog selection

# Generate HIRES scripts

Run 

# Update progress during run. 

cps sync # update logs
run_ckscool.py create-val stat

# Catalog paper

1. Download smemp and smsyn catalogs to data/

$ rsync -av cadence:/data/user/petigura/public_html/smemp/specmatch_results.csv data/specmatch-emp_results.csv
$ rsync -av cadence:/data/user/petigura/public_html/smsyn/specmatch_results.csv data/specmatch-syn_results.csv

2. Run the isoclassify code

```
run_cksgaia.py create-iso-batch 
isoclassify batch direct data/isoclassify-direct.csv -o isoclassify/direct/ > isoclassify-direct.tot
isoclassify batch grid data/isoclassify-grid-parallax-yes.csv  -o isoclassify/grid-parallax-yes/ > isoclassify-grid-parallax-yes.tot
isoclassify batch grid data/isoclassify-grid-parallax-no.csv -o isoclassify/grid-parallax-no/ > isoclassify-grid-parallax-no.tot

# mkdir -p isoclassify/grid-parallax-yes//K00001;isoclassify run grid K00001 --outdir isoclassify/grid-parallax-yes//K00001 --csv data/isoclassify-grid-parallax-yes.csv --dust green18 &> isoclassify/grid-parallax-yes//K00001/output.log

cat isoclassify-direct.tot | parallel 
cat isoclassify-grid-parallax-no.tot | parallel # takes about 20min on Erik's laptop
cat isoclassify-grid-parallax-yes.tot | parallel # takes about 20min on Erik's laptop
```

run_cksgaia.py create-iso-table

3. Generate ReaMatch table

run_cksgaia.py create-reamatch-table 




