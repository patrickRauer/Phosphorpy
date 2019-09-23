import unittest
import pandas as pd

from Phosphorpy.core import structure


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

    def test_getattr(self):
        self.assertTrue((self.table.a == self.data.a).all())
