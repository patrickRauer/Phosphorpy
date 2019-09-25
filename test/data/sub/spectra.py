from astropy import units as u
import numpy as np
import unittest

from Phosphorpy.data.sub.spectra import Spectra, SpectraList


class TestSpectra(unittest.TestCase):

    def setUp(self) -> None:
        self.wavelength = np.float64(np.linspace(360, 900, 900-360))
        self.flux = np.exp(-np.square(self.wavelength-360)/9000)+0.2*np.exp(-np.square(self.wavelength-600)/500)
        self.flux += 0.01*np.random.randn(len(self.flux))
        self.spectra = Spectra(
            wavelength=self.wavelength,
            flux=self.flux
        )

    def test_wavelength_unit(self):
        # check if the default unit is angstrom
        self.assertEqual(
            self.spectra.wavelength_unit,
            u.angstrom
        )

        # if a integer is tried to set as a unit
        with self.assertRaises(ValueError):
            self.spectra.wavelength_unit = 1

        # if a string is tried to set as a unit
        with self.assertRaises(ValueError):
            self.spectra.wavelength_unit = 'srr'

        # if a unit with another dimension is tried to set
        with self.assertRaises(ValueError):
            self.spectra.wavelength_unit = u.deg

        # set nanometer as the new unit
        self.spectra.wavelength_unit = u.nm

        # check if the new units are equivalent
        self.assertEqual(
            self.spectra.wavelength_unit,
            u.nm
        )

    def test_wavelength(self):
        self.assertTrue(
            np.array_equal(
                self.wavelength,
                self.spectra.wavelength
            )
        )

    def test_flux(self):
        self.assertTrue(
            np.array_equal(
                self.flux,
                self.spectra.flux
            )
        )

    def test_cut(self):
        # check, if an error raises, if no limit is set
        with self.assertRaises(ValueError):
            self.spectra.cut(inplace=False)

        # not inplace
        # without units
        spec_cut = self.spectra.cut(min_wavelength=400, inplace=False)
        self.assertGreaterEqual(
            spec_cut.min_wavelength,
            400
        )

        spec_cut = self.spectra.cut(max_wavelength=800, inplace=False)
        self.assertLessEqual(
            spec_cut.max_wavelength,
            800
        )

        spec_cut = self.spectra.cut(min_wavelength=400, max_wavelength=800, inplace=False)
        self.assertGreaterEqual(
            spec_cut.min_wavelength,
            400
        )
        self.assertLessEqual(
            spec_cut.max_wavelength,
            800
        )

        with self.assertRaises(ValueError):
            self.spectra.cut(min_wavelength=400, max_wavelength=800, inplace='hans')

    def test_normalize(self):

        with self.assertRaises(ValueError):
            self.spectra.normalize('summ', inplace=False)

        for kind, func in zip(['sum', 'mean', 'median', 'max'], [np.sum, np.mean, np.median, np.max]):
            norm = func(self.spectra.flux)

            if norm == 0:
                with self.assertRaises(ValueError):
                    self.spectra.normalize(kind, inplace=False)
            else:
                spec_norm = self.spectra.normalize(kind, inplace=False)
                self.assertTrue(
                    np.array_equal(
                        spec_norm.flux,
                        self.spectra.flux/func(self.spectra.flux)
                    )
                )

            if norm == 0:
                with self.assertRaises(ValueError):
                    self.spectra.normalize(kind.upper(), inplace=False)
            else:
                spec_norm = self.spectra.normalize(kind.upper(), inplace=False)
                self.assertTrue(
                    np.array_equal(
                        spec_norm.flux,
                        self.spectra.flux/func(self.spectra.flux)
                    )
                )

    def test_plotting(self):
        self.spectra.plot.spectra()

        self.spectra.plot.spectra(min_wavelength=400)

        self.spectra.plot.spectra(max_wavelength=700)

        self.spectra.plot.spectra(min_wavelength=400, max_wavelength=700)

    def test_fit(self):
        print(self.spectra.fit_line())

    def test_fit_gaus(self):
        print(self.spectra.fit_gauss({'mean': 360, 'stddev': 100}))

    def test_fit_dgaus(self):
        print(self.spectra.fit_double_gauss({'mean': 360, 'stddev': 100},
                                            {'mean': 600, 'stddev': 10}))
