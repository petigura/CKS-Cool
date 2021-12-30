from lmfit import Parameters, Minimizer, report_fit, conf_interval, printfuncs
import ckscool.io
import lmfit
import numpy as np
import pandas as pd

class Fitter(object):
    def __init__(self, mode, ):
        """
        args:
            mode - one of sn-mma - sub-neptune, mass-metallicity-age
                          
        """
        self.mode = mode
        self.nsamples = 1000
        df = ckscool.io.load_table('planets-cuts2')
        df = df[~df.isany]

        self.logage0 = np.log10(5e9)
        self.met0 = 0
        self.logmass0 = 0 
        self.logR_0 = np.log10(2)

        params = Parameters()
        params.add('logR_0', value=self.logR_0)
        params.add('alpha_mass', value=0.00)
        params.add('alpha_met', value=0.0)
        params.add('alpha_age', value=0.0)
        params.add('logmass0', value=self.logmass0)
        params.add('met0', value=self.met0)
        params.add('logage0', value=self.logage0)

        size, fit_method = mode[:2],mode[3:]
        if size=='sn':
            df = df.query('1.7 < gdir_prad < 4.0')
        elif size=='se':
            df = df.query('1.0 < gdir_prad < 1.7')
        else:
            assert False
            
        if fit_method=='m':
            params['alpha_mass'].vary = True 
            params['alpha_met'].vary = False
            params['alpha_age'].vary = False
        elif fit_method=='mm':
            params['alpha_mass'].vary = True 
            params['alpha_met'].vary = True
            params['alpha_age'].vary = False
        elif fit_method=='mma':
            params['alpha_mass'].vary = True 
            params['alpha_met'].vary = True
            params['alpha_age'].vary = True
            df = df.query('cks_steff > 5500')
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
                  + params['alpha_mass'] * (logmass - params['logmass0']) 
                  + params['alpha_met'] * (met - params['met0'])
                  + params['alpha_age'] * (logage - params['logage0']) )
        return _model

    def objective(self, params, df):
        _model = self.model(
            params, np.log10(df.giso_smass), df.cks_smet, df.giso_slogage
        )
        resid = np.log10(df['gdir_prad']) - _model
        return resid 

    def compute_samples(self):
        samples = []
        for i in range(self.nsamples):
            if i %100 ==0:
                print i
            df = self.df.sample(n=self.nplanets,replace=True)
            mini = Minimizer(
                self.objective, self.params0, fcn_args=(df,), nan_policy='omit'
            )
            result = mini.minimize(method='leastsq')
            d = result.params.valuesdict()
            samples.append(d)
         
        samples = pd.DataFrame(samples)
        samples['R0'] = samples.eval('10**logR_0')
        self.samples = samples

