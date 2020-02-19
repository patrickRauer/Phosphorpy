from astropy.table import Table
import numpy as np
import pandas as pd
import os
import shutil

import unittest
try:
    from timeout_decorator import timeout

except ImportError:
    def timeout(a):
        pass

from Phosphorpy.data import data


class TestData(unittest.TestCase):

    def setUp(self) -> None:
        self.image_temp_path = './temp/image/'
        coord = pd.DataFrame(
            {
                'ra': [18.15629, 166.12397, 260.78809],
                'dec': [33.49162, 8.6418, 48.31078]
            }
        )
        Table.from_pandas(coord).write('temp_coordinates.fits', overwrite=True)
        coord.to_csv('temp_coordinates.csv')
        self.ds = data.DataSet.load_coordinates('temp_coordinates.fits')
        # self.ds = data.DataSet()

    def tearDown(self) -> None:
        os.remove('temp_coordinates.fits')
        os.remove('temp_coordinates.csv')
        try:
            shutil.rmtree(self.image_temp_path)
        except:
            pass
        try:
            shutil.rmtree('./temp/')
        except:
            pass

    def test_remove_unmasked_data(self):
        self.ds.remove_unmasked_data()

    def test__get_attribute__(self):
        self.assertTrue(
            (
                self.ds['coordinates'] == self.ds.coordinates
            )
        )

    def test__get_row__(self):
        print('get row')
        print('results', self.ds[1])

    # def test__load_from_vizier__(self):
    #     with self.assertRaises(AttributeError):
    #         self.ds.__load_from_vizier__('hans')

    # def test_load_from_vizier(self):
    #     self.ds.load_from_vizier('all')

    # def test_add_magnitudes(self):
    #     raise NotImplementedError()

    # def test_get_simbad_data(self):
    #     print('simbad')
    #     print(self.ds.get_simbad_data())

    # def test_images(self):
    #     self.ds.images('ps', 1, self.image_temp_path)
    #     list_dir = os.listdir(self.image_temp_path)
    #     self.assertEqual(
    #         len(list_dir), 1
    #     )
    #     shutil.rmtree(self.image_temp_path)
    #
    #     with self.assertRaises(ValueError):
    #         self.ds.images('2MASS', 1, self.image_temp_path)

    # def test_all_images(self):
    #     self.ds.all_images('ps', self.image_temp_path)
    #     list_dir = os.listdir(self.image_temp_path)
    #     self.assertEqual(
    #         len(list_dir), 3
    #     )
    #     shutil.rmtree(self.image_temp_path)

    # def test_light_curves(self):
    #     lc = self.ds.light_curves
    #     print(lc)
    #     print(type(lc))

    # def test_correct_extinction(self):
    #     with self.assertRaises(AttributeError):
    #         self.ds.correct_extinction()
    #     self.ds.load_from_vizier('sdss')
    #     self.ds.correct_extinction()

    def test_reset_masks(self):
        self.ds.reset_masks()

    def test___combine_all__(self):
        print(self.ds.__combine_all__())

    def _write_coordinates_only(self):
        self.ds.write('./temp/ds.zip')

        self.ds.write('./temp/ds2.zip', 'csv')
        self.ds.write('./temp/ds.fits', 'fits')
        self.ds.write('./temp/ds.csv', 'csv')

        with self.assertRaises(ValueError):
            self.ds.write('./temp/ds.txt', 'txt')

    def _write_coordinates_magnitudes1(self):
        self.ds.write('./temp/ds.zip')

        self.ds.write('./temp/ds2.zip', 'csv')
        self.ds.write('./temp/ds.fits', 'fits')
        self.ds.write('./temp/ds.csv', 'csv')

        with self.assertRaises(ValueError):
            self.ds.write('./temp/ds.txt', 'txt')

    def _write_full(self):
        self.ds.write('./temp/ds.zip')

        self.ds.write('./temp/ds2.zip', 'csv')
        self.ds.write('./temp/ds.fits', 'fits')
        self.ds.write('./temp/ds.csv', 'csv')

        with self.assertRaises(ValueError):
            self.ds.write('./temp/ds.txt', 'txt')

    def test_write(self):
        self._write_coordinates_only()
        self.ds.load_from_vizier('sdss')
        self._write_coordinates_magnitudes1()
        self.ds.flux
        self.ds.colors
        self._write_full()

    def test_read_from_file(self):
        with self.assertRaises(FileNotFoundError):
            data.DataSet.read_from_file('test.sdfs')
        self.ds.load_from_vizier('sdss')
        self.ds.flux
        self.ds.colors
        self._write_full()
        ds = data.DataSet.read_from_file('./temp/ds.zip')
        print(ds)

    def test_load_coordinates(self):
        data.DataSet.load_coordinates('temp_coordinates.fits')
        data.DataSet.load_coordinates('temp_coordinates.csv', 'csv')

        with self.assertRaises(ValueError):
            data.DataSet.load_coordinates('temp_coordinates.txt', 'txt')

    # @timeout(1)
    # def test_from_vizier(self):
    #     ds = data.DataSet.from_vizier('SDSS',
    #                                   **{'rmag': '<13', 'spCl': 'STAR'})


class TestData2(unittest.TestCase):

    def setUp(self) -> None:
        self.coord = pd.DataFrame(
            dict(
                ra=[300.706346730, 300.70666426, 300.706019630,
                    300.707633470, 300.70814224, 300.70829226,
                    300.70375521],
                dec=[25.23629731, 25.23574129, 25.237510750,
                     25.237205890, 25.23575684, 25.23468929,
                     25.23600555]
            )
        )
        self.ds = data.DataSet.from_coordinates(self.coord)

    def test_gs2000_ps(self):
        self.ds.load_from_vizier('ps')

        self.assertEqual(
            len(self.coord), len(self.ds.magnitudes['ps'])
        )

    def test_gs2000_all(self):
        self.ds.load_from_vizier('all')
        print(self.ds.magnitudes)
        print(self.ds.magnitudes.data)

        for m in self.ds.magnitudes.data:
            print(m)
            print(m.data)


def test_add_magnitudes():
    coord = pd.DataFrame(
        {
            'ra': [18.15629, 166.12397, 260.78809],
            'dec': [33.49162, 8.6418, 48.31078]
        }
    )
    Table.from_pandas(coord).write('temp_coordinates.fits', overwrite=True)
    ds = data.DataSet.load_coordinates('temp_coordinates.fits')
    ds.load_from_vizier('SDSS')
    ds.load_from_vizier('ps')
    os.remove('temp_coordinates.fits')

    ds.magnitudes.get_survey_data('ps')


def test_load_coordinates():
    data.DataSet.load_coordinates('/Users/patrickr/Documents/test_dash.fits',
                                  'fits', 'RAJ2000', 'DECJ2000')


def test_color_2mass_plot():
    coord = pd.DataFrame(
        {
            'ra': [18.15629, 166.12397, 260.78809],
            'dec': [33.49162, 8.6418, 48.31078]
        }
    )
    Table.from_pandas(coord).write('temp_coordinates.fits', overwrite=True)
    ds = data.DataSet.load_coordinates('temp_coordinates.fits')
    ds.load_from_vizier('SDSS')
    ds.load_from_vizier('ps')
    ds.load_from_vizier('2MASS')
    ds.colors.plot.color_color('2MASS')
