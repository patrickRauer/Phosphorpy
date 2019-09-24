from astropy import units as u


class Spectra:

    _wavelength = None
    _flux = None

    _wavelength_unit = None
    _flux_unit = None

    _lines = None

    def __init__(self):
        pass

    @property
    def wavelength(self):
        return self._wavelength

    @property
    def wavelength_unit(self):
        return self._wavelength_unit

    @wavelength_unit.setter
    def wavelength_unit(self, unit):
        if type(unit) != u.Unit:
            raise ValueError('The new unit must be a astropy unit.')
        else:
            try:
                factor = (self._wavelength_unit.to(unit))
                self._wavelength *= factor
                self._wavelength_unit = unit
            except u.UnitConversionError:
                raise ValueError('Unit must be equivalent to a length.')

    @property
    def flux(self):
        return self._flux

    @property
    def flux_unit(self):
        return self._flux_unit

    @flux_unit.setter
    def flux_unit(self, unit):
        if type(unit) != u.Unit:
            raise ValueError('The new unit must be a astropy unit.')
        else:
            try:
                factor = (self._flux_unit.to(unit))
                self._flux *= factor
                self._flux_unit = unit
            except u.UnitConversionError:
                raise ValueError('Unit must be equivalent to a flux.')
