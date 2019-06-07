sources=(cks1 smsyn smemp)
for source in "${sources[@]}"
do
    isoclassify batch direct data/isoclassify-${source}-direct.csv -o isoclassify/${source}/direct/ > isoclassify-${source}-direct.tot
    isoclassify batch grid data/isoclassify-${source}-grid-parallax-yes.csv  -o isoclassify/${source}/grid-parallax-yes/ > isoclassify-${source}-grid-parallax-yes.tot
    isoclassify batch grid data/isoclassify-${source}-grid-parallax-no.csv -o isoclassify/${source}/grid-parallax-no/ > isoclassify-${source}-grid-parallax-no.tot
done
