from astropy.table import Table
import numpy as np
import pandas as pd
import os
import unittest

from Phosphorpy.data.sub import table


class TestMask(unittest.TestCase):

    def setUp(self) -> None:
        self.mask = table.Mask(10)

    def test_str(self):
        self.assertEqual(str(self.mask),
                         f'Mask No 0: initialization\n')

    def test_add_mask(self):
        with self.assertRaises(ValueError):
            self.mask.add_mask(np.array(10*[True]), 'test')

        self.mask.add_mask(pd.Series(np.array(10*[True])), 'test')

        with self.assertRaises(ValueError):
            self.mask.add_mask(pd.Series(np.array(11*[True])), 'test',
                               combine=False)

        with self.assertRaises(ValueError):
            self.mask.add_mask(pd.Series(np.array(0*[True])), 'test')

    def test_get_latest_mask(self):
        self.assertTrue(
            (self.mask.get_latest_mask() == self.mask.get_mask(-1)).all()
        )

    def test_get_latest_description(self):
        self.assertEqual(
            self.mask.get_latest_description(), self.mask.get_description(-1)
        )

    def test_get_mask(self):
        self.assertTrue(
            (self.mask.get_mask(0) == self.mask._mask[0]).all()
        )

    def test_get_description(self):
        self.assertEqual(
            self.mask.get_description(0),
            self.mask._desc[0]
        )

    def test_get_mask_count(self):
        self.assertEqual(self.mask.get_mask_count(),
                         len(self.mask))

    def test_remove_mask(self):
        length = len(self.mask)
        self.mask.remove_mask(0)
        self.assertEqual(
            len(self.mask),
            length-1
        )

    def test_reset_mask(self):
        self.mask.reset_mask()
        self.assertEqual(len(self.mask), 1)

    def test_mask(self):
        self.assertTrue(
            (self.mask.mask == self.mask.get_latest_mask()).all()
        )

    def test_description(self):
        self.assertTrue(
            (self.mask.description == self.mask.get_latest_description())
        )


class TestDataTable(unittest.TestCase):

    def setUp(self) -> None:
        self.data = pd.DataFrame({'a': [1, 2, 3, 4],
                                  'b': [12., 23., -23, -0.1]},
                                 index=np.arange(4)+1)
        self.table = table.DataTable(self.data.copy())
        # self.table._data = self.data.copy()

    def test_set_mask(self):
        with self.assertRaises(ValueError):
            self.table.set_mask('a')

    def test_stats(self):
        stats = self.table.stats()
        for cols in ['a', 'b']:
            for func, name in zip([np.mean, np.median,
                                   np.std, np.min, np.max],
                                  ['mean', 'median', 'std',
                                   'amin', 'amax']):
                if name == 'std':
                    self.assertEqual(stats[cols][name],
                                     func(self.data[cols].values, ddof=1))
                else:
                    self.assertEqual(stats[cols][name],
                                     func(self.data[cols].values))

    def test_apply(self):
        self.assertTrue(
            (self.table.apply(np.mean) ==
             self.data.apply(np.mean)).all()
        )

    def test_apply_on_ndarray(self):
        self.assertTrue(
            (self.table.apply_on_ndarray(np.mean) ==
             np.mean(self.data.values)).all()
        )

    def test_apply_on_dataframe(self):
        self.assertTrue(
            (self.table.apply_on_dataframe(np.mean) ==
             self.data.apply(np.mean)).all()
        )

    def test_remove_unmask_data(self):
        print()
        mask = self.data['a'] < 3
        print('original mask\n', mask)
        print('default mask\n', self.table.mask.get_latest_mask())
        self.table.mask.add_mask(mask, 'test_mask')
        print('stored mask\n', self.table.mask.get_latest_mask())
        self.table.remove_unmasked_data()
        print(self.table.data)
        print(len(self.table), len(self.table.data), np.sum(self.data['a'].values < 3))
        self.assertEqual(len(self.table), np.sum(self.data['a'] < 3))

    def test_select_nan(self):
        self.table.select_nan('a')

    def test_shape(self):
        shape = self.table.shape
        self.assertEqual(shape[0], len(self.data))
        self.assertEqual(shape[1], len(self.data.columns))

    # def test_head(self):
    #     raise NotImplementedError()

    def test_data(self):
        self.assertTrue(
            (self.table.data == self.data).all().all()
        )

    # def test_plot(self):
    #     raise NotImplementedError()

    def test_q(self):
        self.assertEqual(self.table.q,
                         [0.15, 0.25, 0.75, 0.85])

        self.table.q = 0.2
        self.assertEqual(self.table.q, 0.2)

        self.table.q = [0.15, 0.25, 0.75, 0.85]
        self.assertEqual(self.table.q,
                         [0.15, 0.25, 0.75, 0.85])

        with self.assertRaises(ValueError):
            self.table.q = [0.15, 0.25, 0.75, 0.85, -1]

        with self.assertRaises(ValueError):
            self.table.q = [0.15, 0.25, 0.75, 0.85, 2]

    def test_mask(self):
        self.assertTrue(
            (self.table.mask == self.table._mask)
        )

    def test_to_astropy_table(self):
        tab = self.table.to_astropy_table()
        self.assertTrue(
            (Table.from_pandas(self.data) == tab).all()
        )

    def test_write(self):

        for f in ['csv', 'latex', 'fits', 'sql', 'parquet']:
            path = f'./temp.{f}'
            self.table.write(path, f)
            self.assertTrue(os.path.exists(path))
            os.remove(path)

        with self.assertRaises(ValueError):
            self.table.write('./temp.txt', 'txt')

    def test_str(self):
        self.assertEqual(
            str(self.table), str(self.data)
        )

    def test_getitem(self):
        self.assertTrue(
            (self.table['a'] == self.data['a']).all()
        )


if __name__ == '__main__':
    unittest.main()
