"""
Module to specify CKS-Cool sample

- exclude stars with known dillution from Furlan paper? How many more would we get?
- Borrow code from CKS-cuts.
"""



cks = cksphys.io.load_table('cksphys-merged',cache=1)
stellar17 = cksphys.io.load_table('stellar17',cache=1)
nea = cksphys.io.load_table('nea-cum',cache=1)
nea = nea[nea.koi_disposition.str.contains('C')]
df = pd.merge(nea,stellar17,on='id_kic')
cks['incks'] = 1
df = pd.merge(df,cks['id_kic incks'.split()],how='left')
df = df.fillna(0)

# Add in Kraus AO
import astropy.io.ascii
from astropy.table import Table
t = Table.read('../NASA-2018A/kraus16/table2.dat',readme='../NASA-2018A/kraus16/ReadMe',format='ascii.cds')
kraus = t.to_pandas()
kraus['inkraus16'] = 1
kraus = kraus.rename(columns={"KOI":"id_koi"})
df = pd.merge(df,kraus['id_koi inkraus16'.split()],how='left')
df = df.fillna(0)

tcand = df.query('kic_steff < 4900 and kic_kepmag < 16')
tstar = tcand.groupby('id_kic',as_index=False).first()

ckscand = cks
cksstar = ckscand.groupby('id_kic',as_index=False).first()

len(tcand)
len(tstar)
