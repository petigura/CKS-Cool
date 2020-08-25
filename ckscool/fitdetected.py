from lmfit import Parameters, Minimizer, report_fit, conf_interval, printfuncs
import ckscool.io
import lmfit
import numpy as np
import pandas as pd

class Fitter(object):
    def __init__(self, mode, ):
        self.nsamples = 1000
        df = ckscool.io.load_table('planets-cuts2')
        df = df[~df.isany]
        params = Parameters()
        params.add('logR_0', value=np.log10(2))
        params.add('alpha_mass', value=0.00)
        params.add('alpha_met', value=0.0)
        params.add('alpha_age', value=0.0)
        params.add('logmass0', value=0)
        params.add('met0', value=0)
        params.add('logage0', value=np.log10(5e9))

        size, fit_method = mode[:2],mode[3:]
        if size=='sn':
            df = df.query('1.7 < gdir_prad < 4.0')
        elif size=='se':
            df = df.query('1.0 < gdir_prad < 1.7')
        else:
            assert False
            
        if fit_method=='mass':
            params['alpha_mass'].vary = True 
            params['alpha_met'].vary = False
            params['alpha_age'].vary = False
        elif fit_method=='mass-met':
            params['alpha_mass'].vary = True 
            params['alpha_met'].vary = True
            params['alpha_age'].vary = False
        elif fit_method=='mass-met-age':
            params['alpha_mass'].vary = True 
            params['alpha_met'].vary = True
            params['alpha_age'].vary = True
        else:
            assert False

        self.params0 = params
        self.df = df
        self.nplanets = len(df)

    def model(self, params, logmass, met, logage):
        """
        params: lmfit parameter object
        logmass: pandas dataframe log10 of stellar mass
        met (pd.DataFrame): metallicity
        logage (pd.DataFrame): log of planet mass
        """

        _model = (params['logR_0'] 
                  + params['alpha_mass']*(logmass-params['logmass0']) 
                  + params['alpha_met']*(met-params['met0'])
                  + params['alpha_age']*(logage-params['logage0']) )
        return _model

    def objective(self, params, df):
        _model = self.model(
            params, np.log10(df.giso_smass), df.cks_smet,df.giso_slogage
        )
        resid = np.log10(df['gdir_prad']) - _model
        return resid 

    def compute_samples(self):
        samples = []
        for i in range(self.nsamples):
            if i %100 ==0:
                print i
            df = self.df.sample(n=self.nplanets,replace=True)
            mini = Minimizer(self.objective, self.params0, fcn_args=(df,), nan_policy='omit')
            result = mini.minimize(method='leastsq')
            d = result.params.valuesdict()
            samples.append(d)
            
        samples = pd.DataFrame(samples)
        self.samples = samples


'''
df = df.query('koi_period < 300')


logmassi = linspace(log10(0.5),log10(1.4),100)
meti = zeros_like(logmassi)
logagei = zeros_like(logmassi) + 9 
for i,row in dL.sample(100).iterrows():
    params = Parameters()
    for k in 'logR_0 alpha_mass alpha_met alpha_age logmass0 met0 logage0'.split():
        params.add(k,value=row[k])
    _model = model(params, logmassi, meti, logagei)
    loglog()
    plot(10**logmassi,10**_model,color='r')
plot(cut.giso_smass,cut.gdir_prad,'.') 
ylabel('Stellar Mass')
xlabel('Metallicity')
tight_layout()

clf()
df = ckscool.io.load_table('planets-cuts2')
plot(df.giso_smass,df.giso_slogage,'.')
df = df.query('~isany')
plot(df.giso_smass,df.giso_slogage,'.')
plot(df.giso_smass,df.giso_slogage,'.')

%pylab osx
import ckscool.io
import lmfit
df = ckscool.io.load_table('planets-cuts2')
df = df[~df.isany]
#df = df.query('0.6 < giso_smass < 1.2 and -0.2 < cks_smet < 0.2' )
df = df.query('koi_period < 300')
df = df.query('cks_steff  >5500')

#cut0 = df.query('1.7 < gdir_prad < 4.0')
cut0 = df.query('1.0 < gdir_prad < 1.7')
#fracdev = 10**np.log10(cut.gdir_prad).std() - 1
err = cut.gdir_prad.std()

from lmfit import Parameters, Minimizer, report_fit, conf_interval, printfuncs
params = Parameters()
params.add('logR_0', value=log10(2),vary=True)
params.add('alpha_mass', value=0.00,vary=True)
params.add('alpha_met', value=0.0,vary=True)
params.add('alpha_age', value=0.0,vary=True)
params.add('logmass0', value=0,vary=False)
params.add('met0', value=0,vary=False)
params.add('logage0', value=9,vary=False)

def model(params, logmass, met, logage):
    _model = (params['logR_0'] 
              + params['alpha_mass']*(logmass - params['logmass0']) 
              + params['alpha_met']*(met - params['met0'])
              + params['alpha_age']*(logage - params['logage0']) )
    return _model

cut = cut0.copy()
def objective(params):
    _model = model(params, np.log10(cut['giso_smass']),cut['cks_smet'],cut['giso_slogage'])
    resid = log10(cut['gdir_prad']) - _model
    return resid 

mini = Minimizer(objective, params, nan_policy='omit')
result = mini.minimize(method='Nelder')
result = mini.minimize(method='leastsq',params=params)
print(report_fit(result))

dL = []
for i in range(1000):
    d = result.params.valuesdict()
    if i %100 ==0:
        print i
    cut = cut0.sample(n=len(cut),replace=True)
    def objective(params):
        _model = model(params, np.log10(cut['giso_smass']),cut['cks_smet'],cut['giso_slogage'])
        resid = log10(cut['gdir_prad']) - _model
        return resid 

    mini = Minimizer(objective, params, nan_policy='omit')
    result = mini.minimize(method='leastsq',params=result.params)
    dL.append(d)

dL = pd.DataFrame(dL)

clf()
df = ckscool.io.load_table('planets-cuts2')
plot(df.giso_smass,df.giso_slogage,'.')
df = df.query('~isany')
plot(df.giso_smass,df.giso_slogage,'.')
plot(df.giso_smass,df.giso_slogage,'.')

clf()
plot(df.giso_smass, df.giso_slogage,'.')
'''
