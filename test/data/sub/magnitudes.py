import numpy as np
import pandas as pd
import os
import unittest
import warnings
from urllib.error import URLError

from Phosphorpy.data.sub import magnitudes


class TestMagnitudes(unittest.TestCase):

    def test_mag_cols_only(self):
        cols = ['gmag', 'Kmag', 'e_gmag', 'e_Kmag']
        self.assertEqual(magnitudes.mag_cols_only(cols),
                         ['gmag', 'Kmag'])


class TestSurveyData(unittest.TestCase):

    def setUp(self) -> None:
        try:
            self.survey = magnitudes.SurveyData('2mass',
                                                ['J', 'H', 'K'])
        except URLError:
            self.skipTest('Could not test survey properties because the SVO web-service is down.')

    def test__add_survey_properties__(self):
        try:
            self.survey.__add_survey_properties__(
                'SDSS', ['u', 'g']
            )
            with self.assertRaises(ValueError):
                self.survey.__add_survey_properties__(
                    'SDSS', ['u', 'v']
                )
            with self.assertRaises(KeyError):
                self.survey.__add_survey_properties__(
                    'HST', ['u', 'v']
                )
        except URLError:
            warnings.warn('Could not test survey properties because the SVO web-service is down.')

    def test_getitem(self):
        self.assertEqual(
            self.survey['2mass'],
            self.survey._survey_cols['2mass']
        )
        self.assertEqual(
            self.survey['2MASS'],
            self.survey._survey_cols['2mass']
        )

    def test_all_magnitudes(self):
        self.assertEqual(
            self.survey.all_magnitudes(),
            ['J', 'H', 'K']
        )

    def test_get_survey_magnitude(self):
        self.assertEqual(
            self.survey.get_survey_magnitude('2mass'),
            ['J', 'H', 'K']
        )
        self.assertEqual(
            self.survey.get_survey_magnitude('2MASS'),
            ['J', 'H', 'K']
        )

    def test_get_survey_error_names(self):
        self.assertEqual(
            self.survey.get_survey_error_names('2mass'),
            ['e_J', 'e_H', 'e_K']
        )
        self.assertEqual(
            self.survey.get_survey_error_names('2MASS'),
            ['e_J', 'e_H', 'e_K']
        )

    def test_all_error_names(self):
        self.assertEqual(
            self.survey.all_error_names(),
            ['e_J', 'e_H', 'e_K']
        )

    def test_get_surveys(self):
        self.assertEqual(
            self.survey.get_surveys(),
            ['2mass']
        )
        self.assertEqual(
            self.survey.get_surveys(),
            ['2mass']
        )

    def test_add_survey(self):
        self.survey.add_survey(
            'wise', ['W1', 'W2', 'W3', 'W4']
        )
        self.survey.add_survey(
            'WISE', ['W1', 'W2', 'W3', 'W4']
        )

    def test_effective_wavelength(self):
        self.assertEqual(
            self.survey.effective_wavelength(
                '2mass', 'K'
            ),
            21521.64
        )
        self.assertEqual(
            self.survey.effective_wavelength(
                '2MASS', 'K'
            ),
            21521.64
        )

    def test_get_all_wavelengths(self):
        self.assertTrue(
            (self.survey.get_all_wavelengths() ==
             np.array([12285.38, 16386.10, 21521.64])).all()
        )

    def test_get_survey_wavelengths(self):
        self.assertTrue(
            (self.survey.get_survey_wavelengths('2mass') ==
             np.array([12285.38, 16386.10, 21521.64])).all()
        )
        self.assertTrue(
            (self.survey.get_survey_wavelengths('2MASS') ==
             np.array([12285.38, 16386.10, 21521.64])).all()
        )

    def test_to_dataframe(self):
        self.assertTrue((
            self.survey.to_dataframe().columns.values ==
            np.array(['survey', 'band', 'lambda_mean', 'lambda_cen', 'lambda_eff',
                      'lambda_peak', 'lambda_pivot', 'lambda_phot', 'lambda_min',
                      'lambda_max', 'w_eff', 'fwhm', 'af_av', 'vega_ergs', 'vega_jy',
                      'ab_ergs', 'ab_jy', 'st_ergs', 'st_jy', 'vega'])).all()
        )

    def test_write(self):
        self.survey.write('temp.ini')
        self.assertTrue(
            os.path.exists('temp.ini')
        )
        os.remove('temp.ini')

        self.survey.write('temp.sql', 'sql')
        self.assertTrue(
            os.path.exists('temp.sql')
        )
        os.remove('temp.sql')


class TestMagnitudeTable(unittest.TestCase):

    def setUp(self) -> None:
        self.data = pd.DataFrame(
            {'umag': [0, 0, 0],
             'gmag': [1, 2, 3],
             'rmag': [1.1, 2.2, 3.3],
             'imag': [-1, -2, -3],
             'zmag': [-1.1, -2.2, -3.3],
             'e_umag': [0, 0, 0],
             'e_gmag': [1, 2, 3],
             'e_rmag': [1.1, 2.2, 3.3],
             'e_imag': [-1, -2, -3],
             'e_zmag': [-1.1, -2.2, -3.3]
             }
        )
        self.mags = magnitudes.MagnitudeTable(

        )
        self.mags.add_survey_mags(self.data.copy(),
                                  'SDSS')

    def test_survey(self):
        print(self.mags.survey)

    def test_create_mask(self):
        self.mags.create_mask()

        with self.assertRaises(RuntimeError):
            self.mags.create_mask()

        self.mags.create_mask(overwrite=True)

    def test_get_names(self):
        self.assertEqual(
            self.mags.get_names(),
            ['SDSS']
        )

    def test_get_flux(self):
        flux = self.mags.get_flux()

    def test_get_colors(self):
        colors = self.mags.get_colors()

    def test_has_full_photometry(self):
        self.mags.has_full_photometry('SDSS')

    def test_set_limit(self):
        self.mags.set_limit('u', 'SDSS', 1, 2)

    def test_get_survey_data(self):
        sd = self.mags.get_survey_data('SDSS')

        sd = self.mags.get_survey_data('sdss')

    def test_get_columns(self):
        cols = ['gmag', 'rmag']
        d = self.mags.get_columns(cols)
        self.assertEqual(cols, list(d.columns.values))

        self.assertIsNone(self.mags.get_columns(['Kmag']))

        self.assertEqual(['umag'], list(self.mags.get_columns('umag')))

        with self.assertRaises(AttributeError):
            self.mags.get_columns(1)


if __name__ == '__main__':
    unittest.main()