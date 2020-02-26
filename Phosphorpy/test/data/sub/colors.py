import numpy as np
import pandas as pd

import unittest

from Phosphorpy.data.sub import colors


class TestColors(unittest.TestCase):

    def setUp(self) -> None:
        self.data = pd.DataFrame(
            {
                'g-r': np.random.random(10),
                'r-i': np.random.random(10)
            }
        )
        self.color = colors.ColorTab(self.data.copy(), 'sdss')
        self.colors = colors.Colors([self.color])

    def test_add_colors(self):
        data = pd.DataFrame(
            {
                'g-r': np.random.random(10),
                'r-i': np.random.random(10)
            }
        )
        self.colors.add_colors(data.copy(), 'panspanstarrstar')

    def test__get_mask_data__(self):
        self.colors.__get_mask_data__(
            'g-r', 0, 1, True
        )

    def test_set_limit(self):
        self.colors.set_limit('g-r', 0, 2)

        with self.assertRaises(ValueError):
            self.colors.set_limit('g-r', 2, 0)

    # def test_get_columns(self):
    #     cols = self.colors.get_columns('g-r')
    #     self.assertTrue(
    #         (
    #             cols == self.data[['g-r']]
    #         ).all().all()
    #     )

    def test_get_column(self):
        col = self.colors.get_column('g-r')
        self.assertTrue(
            (
                col == self.data['g-r']
            ).all()
        )

    def test_outlier_detection(self):
        self.colors.outlier_detection('sdss')

    def test_str(self):
        str(self.colors)

    def test_survey_colors(self):
        self.colors.survey_colors


if __name__ == '__main__':
    unittest.main()