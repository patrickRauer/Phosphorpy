import numpy as np
import pandas as pd
import unittest

from Phosphorpy.data.sub.tables.flux import Flux


class TestFlux(unittest.TestCase):

    def setUp(self) -> None:
        data = pd.DataFrame(
            {'umag': np.random.rand(10),
             'e_umag': np.random.rand(10),
             'gmag': -np.random.rand(10),
             'e_gmag': np.random.rand(10)}
        )
        self.flux = Flux(data.copy())
        self.data = data

    def test_get_fluxes(self):
        d = self.flux.get_fluxes()
        self.assertEqual(['umag', 'gmag'],
                         list(d.columns.values))
        self.assertEqual(len(d), len(self.data))

    def test_get_errors(self):
        d = self.flux.get_errors()
        self.assertEqual(['e_umag', 'e_gmag'],
                         list(d.columns.values))
        self.assertEqual(len(d), len(self.data))

    def test_get_flux(self):
        for index in range(len(self.data)):
            f = self.flux.get_flux(index)
            self.assertTrue(
                (self.data[['umag', 'gmag']].iloc[index] ==
                 f).all()
            )

    def test_get_error(self):
        for index in range(len(self.data)):
            f = self.flux.get_error(index)
            self.assertTrue(
                (self.data[['e_umag', 'e_gmag']].iloc[index] ==
                 f).all()
            )

    def test_get_wavelengths(self):
        self.flux.get_wavelengths()
