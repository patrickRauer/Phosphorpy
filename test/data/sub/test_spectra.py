from astropy import units as u
import numpy as np
import os
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
        self.spectra.plot.spectra()

    def test_fit_dgaus(self):
        print(self.spectra.fit_double_gauss({'mean': 360, 'stddev': 100},
                                            {'mean': 600, 'stddev': 10}))
        print(self.spectra.fit_gauss({'mean': 360, 'stddev': 100}))
        print(self.spectra.fit_gauss({'mean': 600, 'stddev': 10}))
        self.spectra.plot.spectra()

    def test_str(self):
        print(self.spectra)

    def test_write(self):
        self.spectra.write('./test_spec.fits')
        self.spectra.write('./test_spec.csv', data_format='csv')

        with self.assertRaises(ValueError):
            self.spectra.write('./test_spec.hans', data_format='hans')

        os.remove('./test_spec.fits')
        os.remove('./test_spec.csv')

    def test_read(self):
        self.spectra.write('./test_spec.fits')
        Spectra.read('./test_spec.fits')
        self.spectra.write('./test_spec.csv', 'csv')
        Spectra.read('./test_spec.csv', 'csv')

        with self.assertRaises(ValueError):
            Spectra.read('./test_spec.csv', 'hans')

        os.remove('./test_spec.fits')
        os.remove('./test_spec.csv')

    def test_copy(self):
        spec = self.spectra.copy()


class TestSpectraList(unittest.TestCase):

    def setUp(self) -> None:
        self.spec_list = SpectraList()

        wavelength = np.float64(np.linspace(360, 900, 900-360))
        flux = np.exp(-np.square(wavelength-360)/9000)+0.2*np.exp(-np.square(wavelength-600)/500)
        flux += 0.01*np.random.randn(len(flux))
        spectra = Spectra(
            wavelength=wavelength,
            flux=flux
        )
        self.spectra = spectra

    def test_append(self):
        self.spec_list.append(self.spectra)

        self.assertEqual(len(self.spec_list), 1)

        self.spec_list.append(self.spectra, spec_id=2)

        self.assertEqual(len(self.spec_list), 2)

    def test_print(self):
        print(self.spec_list)

    def test_get_item(self):
        self.spec_list.append(self.spectra)
        self.spec_list.append(self.spectra)
        spec, id = self.spec_list[0]

        specs, ids = self.spec_list[[0, 1]]

    def test_merge(self):
        self.spec_list.append(self.spectra)
        self.spec_list.append(self.spectra)
        self.spec_list.merge(self.spec_list)

        with self.assertRaises(ValueError):
            self.spec_list.merge(self.spectra)

    def test_get_by_id(self):
        self.spec_list.append(self.spectra)
        self.spec_list.append(self.spectra)
        self.spec_list.merge(self.spec_list)

        sl = self.spec_list.get_by_id(0)

        sl = self.spec_list.get_by_id(np.array([0, 1]))

    def test_write(self):
        self.spec_list.append(self.spectra)
        self.spec_list.write('./test/')

        with self.assertRaises(ValueError):
            self.spec_list.write('./test/None_0.fits')

    def test_read(self):
        self.spec_list.append(self.spectra)
        self.spec_list.write('./test/')

        SpectraList.read('./test/')

    def test_spectra_list_init(self):
        wavelength = np.float64(np.linspace(360, 900, 900 - 360))
        flux = np.exp(-np.square(wavelength - 360) / 9000) + 0.2 * np.exp(-np.square(wavelength - 600) / 500)
        flux += 0.01 * np.random.randn(len(flux))
        spectra = Spectra(
            wavelength=wavelength,
            flux=flux
        )
        sl = SpectraList()

        sl = SpectraList(spectra)

        sl = SpectraList(spectra, 0)

        sl = SpectraList((spectra, spectra))
        sl = SpectraList((spectra, spectra), [1, 2])

        with self.assertRaises(ValueError):
            sl = SpectraList((spectra, spectra), [0, 1, 2])


if __name__ == '__main__':
    unittest.main()
