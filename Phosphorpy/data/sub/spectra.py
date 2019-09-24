from astropy import units as u


class Spectra:

    _wavelength = None
    _flux = None

    _wavelength_unit = None
    _flux_unit = None

    _lines = None

    def __init__(self, wavelength, flux, wavelength_unit=None, flux_unit=None):
        self._wavelength = wavelength
        self._flux = flux

        self.wavelength_unit = wavelength_unit
        self.flux_unit = flux_unit

    @property
    def wavelength(self):
        return self._wavelength

    @property
    def wavelength_unit(self):
        return self._wavelength_unit

    @wavelength_unit.setter
    def wavelength_unit(self, unit):
        if unit is None:
            self._wavelength_unit = u.angstrom
        elif type(unit) != u.Unit:
            raise ValueError('The new unit must be a astropy unit.')
        elif self._wavelength_unit is None:
            self._wavelength_unit = unit
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
        if unit is None:
            self._flux_unit = u.erg
        elif type(unit) != u.Unit:
            raise ValueError('The new unit must be a astropy unit.')
        elif self._flux_unit is None:
            self._flux_unit = unit
        else:
            try:
                factor = (self._flux_unit.to(unit))
                self._flux *= factor
                self._flux_unit = unit
            except u.UnitConversionError:
                raise ValueError('Unit must be equivalent to a flux.')
