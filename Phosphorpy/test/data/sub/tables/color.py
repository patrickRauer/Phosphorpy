import numpy as np
import pandas as pd

import unittest

from Phosphorpy.data.sub.tables.color import Color


class TestColor(unittest.TestCase):

    def setUp(self) -> None:
        data = pd.DataFrame(
            {
                'g-r': np.random.random(10),
                'r-i': np.random.random(10)
            }
        )
        self.color = Color(data.copy(), 'color')
        self.data = data

    # def test_get_columns(self):
        #     #     self.assertTrue(
        #     #         (self.color.get_columns(['g-r']) == self.data[['g-r']]).all().all()
        #     #     )
        #     #     self.assertTrue(
        #     #         (self.color.get_columns('g-r') == self.data[['g-r']]).all().all()
        #     #     )
        #     #     self.assertTrue(
        #     #         (self.color.get_columns(['g']) == self.data[['g-r']]).all().all()
        #     #     )
        #     #     self.assertTrue(
        #     #         (self.color.get_columns('g') == self.data[['g-r']]).all().all()
        #     #     )
        #     #     self.assertTrue(
        #     #         (self.color.get_columns(['r']) == self.data).all().all()
        #     #     )

        with self.assertRaises(ValueError):
            self.color.get_columns('j-k')

        with self.assertRaises(AttributeError):
            self.color.get_columns(1)

        with self.assertRaises(ValueError):
            self.color.get_columns([1])

    def test_set_limit(self):
        with self.assertRaises(ValueError):
            self.color.set_limit(['g-r'])
        with self.assertRaises(ValueError):
            self.color.set_limit(['g-r', 'r-i', 'i-z', 'z-y'])

        with self.assertRaises(ValueError):
            self.color.set_limit(1)

        self.color.set_limit('g-r')


