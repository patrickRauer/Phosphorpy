

class Spectra:

    _wavelength = None
    _flux = None

    _lines = None

    def __init__(self):
        pass

    @property
    def wavelength(self):
        return self._wavelength

    @property
    def flux(self):
        return self._flux
