import unittest
import pandas as pd
import os

from Phosphorpy.core import structure


class TestMask(unittest.TestCase):

    def setUp(self) -> None:
        self.mask = structure.Mask(2)

    def test_initialization(self):
        mask = structure.Mask(2)
        mask = structure.Mask(2, [1, 2])

        with self.assertRaises(AttributeError):
            structure.Mask(1, [1, 2])

    def test_str(self):
        print(self.mask)

    def test_add_mask(self):
        ds = pd.Series([True, False])
        with self.assertRaises(ValueError):
            self.mask.add_mask(1)

        with self.assertRaises(ValueError):
            self.mask.add_mask(pd.Series())

        self.mask.add_mask(ds)

    def test_get_latest_mask(self):
        self.assertTrue(
            (self.mask.get_latest_mask() == self.mask._mask[-1]).all()
        )

    def test_get_latest_description(self):
        self.assertTrue(
            self.mask.get_latest_description() == self.mask.description
        )

    def test_get_mask(self):
        self.assertTrue(
            (self.mask.get_mask(0) == self.mask._mask[0]).all()
        )

    def test_get_description(self):
        self.assertTrue(
            (self.mask.get_description(0) == self.mask._desc[0])
        )

    def test_remove_mask(self):
        self.mask.remove_mask(0)
        self.assertEqual(len(self.mask), 0)

    def test_reset_mask(self):
        ds = pd.Series([True, False])
        self.mask.add_mask(ds)
        self.mask.reset_mask()
        self.assertEqual(len(self.mask), 1)


class TestHead(unittest.TestCase):

    def setUp(self) -> None:
        self.head = structure.Head()

    def test_name(self):
        self.assertEqual(self.head.name, '')
        self.head.name = 'test'
        self.assertEqual(self.head.name, 'test')

        with self.assertRaises(ValueError):
            self.head.name = None

        with self.assertRaises(ValueError):
            self.head.name = 1


class TestTable(unittest.TestCase):

    def setUp(self) -> None:
        self.data = pd.DataFrame({'a': [1, 2, 3], 'b': [9, 8, 7]})
        self.table = structure.Table(self.data, 'test', None)

    def test_length(self):
        self.assertEqual(len(self.table), len(self.data))

    def test_getitem(self):
        for c in self.data.columns:
            self.assertTrue((self.table[c] == self.data[c]).all())

    def test_write(self):
        path_fits = self.table.write('./test.fits', 'fits')
        path_csv = self.table.write('./test.csv', 'csv')
        path_parquet = self.table.write('./test.parquet', 'parquet')

        with self.assertRaises(ValueError):
            self.table.write('./test.hans', 'hans')

        os.remove(path_fits)
        os.remove(path_csv)
        os.remove(path_parquet)

    def test_rename(self):
        self.table.rename({'a': 'a_new', 'b': 'b_new'}, 'columns')

    def test_select_columns(self):
        self.table.select_columns(['a'])

    def test_merge(self):
        self.table.merge(self.data)


if __name__ == '__main__':
    unittest.main()