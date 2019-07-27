import numpy as np
import xarray as xr 


class Grid(object):
    """
    Grid

    This object creates grids used to compute occurrence and completeness
    """

    def __init__(self, bins_dict, spacing_dict):
        """Lay down grid for computing occurrence
        
        Args:
            bins_dict (dict): Dictionary of lists defining bins. e.g.:
                {'per': [5,10,20], 'prad': [1,2,4], 'smet': [-0.4,0.0,0.4]}
            spacing_dict (dict): Specify linear or log spacing e.g.:
                {'per': 'log', 'prad': 'log', 'smet': 'linear'}
        
        """

        assert bins_dict.keys()==spacing_dict.keys()
        self.bins = bins_dict
        self.spacing = spacing_dict

        # For every bin, record the lower, upper, and central value
        self.bins1 = {}
        self.bins2 = {}
        self.binsc = {}
        self.nbins = {}

        for key in self.bins.keys():
            self._add_bin_edges(key)
            self.nbins[key] = len(self.bins1[key])


        ds = xr.Dataset(coords=self.binsc)
        for sbin in self.bins.keys():
            ds.coords[sbin] = self.binsc[sbin]

        for key in self.bins.keys():
            kw = dict(coords=[self.binsc[key]],dims = [key])
            ds[key+'c'] = xr.DataArray(self.binsc[key],**kw)
            ds[key+'1'] = xr.DataArray(self.bins1[key],**kw)
            ds[key+'2'] = xr.DataArray(self.bins2[key],**kw)

        # This line makes all data variables 2D
        ds = ds.to_dataframe().to_xarray()
        self.ds = ds

    def _add_bin_edges(self, key):
        """Computes bin edges and centers for the specified bins"""
        bins = self.bins[key]
        spacing = self.spacing[key]

        bins1 = []
        bins2 = []
        binsc = []
        for i in range(len(bins)-1):
            bin1 = bins[i]
            bin2 = bins[i+1]
            if spacing=='linear':
                binc = 0.5 * (bin1 + bin2)
            elif spacing=='log':
                binc = np.sqrt(bin1 * bin2)
            else:
                False, "Invalid spacing"
            bins1+=[bin1]
            bins2+=[bin2]
            binsc+=[binc]
        
        self.bins1[key] = bins1
        self.bins2[key] = bins2
        self.binsc[key] = binsc
