from astroquery.irsa_dust import IrsaDust
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import pandas as pd
import warnings


def get_extinctions(ra, dec, wavelengths=None, filter_names=None):
    """
    Fall back method if the extinction package is not available.
    For performance reasons, the extinction in the closest filter is used.

    :param ra: The RA component of the coordinates
    :type ra: Union
    :param dec: The Dec component of the coordinates
    :type dec: Union
    :param wavelengths: The wavelengths of the filter
    :type wavelengths: Union
    :param filter_names: The name of the filters
    :type filter_names: Union
    :return: The extinctions of the sources in the different bands
    :rtype: pandas.DataFrame
    """
    coords = SkyCoord(ra*u.deg, dec*u.deg)
    warnings.warn('Fall back to direct IRSA dust map queries.\nSlow!!!')

    wavelengths /= 1000
    out = []
    for c in coords:
        extinctions = IrsaDust.get_extinction_table(c)
        r = []
        for w, names in zip(wavelengths, filter_names):
            z = np.argmin(np.abs(extinctions['LamEff']-w))
            r.append(extinctions['A_SandF'][z])
        out.append(r)
    return pd.DataFrame(data=np.array(out), columns=filter_names)
