import unittest
import numpy as np
from Phosphorpy.core import functions


class TestFunctions(unittest.TestCase):
    
    def test_smooth(self):
        a = np.array([0., 0., 1., 0., 0.])
        b = np.array([0., 0.25, 0.5, 0.25, 0.])

        a_smoothed = functions.smooth(a, 1)
        a_2smoothed = functions.smooth(a, 2)
        
        self.assertTrue(np.all(b == a_smoothed))
        self.assertTrue(np.all(a_2smoothed == functions.smooth(a_smoothed, 1)))
        
    def test_smooth2d(self):
        a = np.zeros([5, 5])
        a[2, 2] = 1
        
        b = np.zeros((5, 5))
        b[2, 2] = 0.5
        b[2, 1] = 0.125
        b[2, 3] = 0.125
        b[1, 2] = 0.125
        b[3, 2] = 0.125
        
        a_smooth = functions.smooth2d(a, 1)
        
        self.assertTrue(np.all(b == a_smooth))
