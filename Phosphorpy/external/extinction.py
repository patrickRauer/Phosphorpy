from astroquery.irsa_dust import IrsaDust
from astropy.coordinates import SkyCoord
from astropy.table import vstack
from astropy import units as u
import numpy as np
import pandas as pd
import warnings
import numbers


def _to_array(x):
    if type(x) != np.ndarray:
        if isinstance(x, numbers.Number):
            x = np.array([x])
        else:
            x = np.array(x)
    return x


def get_extinctions(ra, dec, wavelengths=None, filter_names=None):
    """
    Fall back method if the extinction package is not available.
    For performance reasons, the extinction in the closest filter is used.

    :param ra: The RA component of the coordinates
    :type ra: Union
    :param dec: The Dec component of the coordinates
    :type dec: Union
    :param wavelengths: The wavelengths of the filter in um or as a astropy Quantity with the dimension length
    :type wavelengths: Union, Quantity
    :param filter_names: The name of the filters
    :type filter_names: Union
    :return: The extinctions of the sources in the different bands
    :rtype: pandas.DataFrame
    """
    warnings.warn('Fall back to direct IRSA dust map queries.\nSlow!!!')
    ra = _to_array(ra)
    dec = _to_array(dec)

    coords = SkyCoord(ra*u.deg, dec*u.deg)

    if wavelengths is None:
        out = []
        for c in coords:
            extinctions = IrsaDust.get_extinction_table(c)
            extinctions['ra'] = c.ra.degree
            extinctions['dec'] = c.dec.degree
            out.append(extinctions)
        return vstack(out).to_pandas()
    else:
        if filter_names is None:
            raise ValueError('If wavelengths are given, filter_names must be given, too.')
        if type(wavelengths) == u.Quantity:
            wavelengths = wavelengths.to(u.um).value

        if isinstance(wavelengths, numbers.Number):
            wavelengths = [wavelengths]
            filter_names = [filter_names]
        out = []
        for c in coords:
            extinctions = IrsaDust.get_extinction_table(c)
            r = []
            for w, names in zip(wavelengths, filter_names):
                z = np.argmin(np.abs(extinctions['LamEff']-w))
                r.append(extinctions['A_SandF'][z])
            out.append(r)
        return pd.DataFrame(data=np.array(out), columns=filter_names)
