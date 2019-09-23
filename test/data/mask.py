import numpy as np
import unittest

from Phosphorpy.data import mask


class TestMaskLog(unittest.TestCase):

    def setUp(self) -> None:
        self.test_mask = np.array([True, False, True, False])
        self.test_label = 'test1'
        self.mask_log = mask.MaskLog(self.test_mask,
                                     self.test_label)

    def test_get_mask(self):
        self.assertTrue((self.mask_log.mask == self.test_mask).all())

    def test_set_mask(self):
        with self.assertRaises(AttributeError):
            self.mask_log.mask = 1

    def test_get_label(self):
        self.assertEqual(self.mask_log.label,
                         self.test_label)

    def test_set_label(self):

        test_label2 = 'test2'
        self.mask_log.label = test_label2
        self.assertEqual(self.mask_log.label,
                         test_label2)


class TestMask(unittest.TestCase):

    def setUp(self) -> None:
        self.test_mask = np.array([True, False, True, False])
        self.test_label = 'test1'
        self.mask_log = mask.MaskLog(self.test_mask,
                                     self.test_label)

        self.test_mask = np.array([True, False, True, False])
        self.test_label = 'test1'
        self.mask_log = mask.MaskLog(self.test_mask,
                                     self.test_label)
        self.mask = mask.Mask()

    def test_set_mask(self):
        self.mask.mask = self.test_mask

    def test_setitem(self):
        self.mask['new mask'] = self.test_mask

    def test_get_mask(self):
        with self.assertRaises(ValueError):
            self.mask.mask

        self.mask.mask = self.test_mask
        self.assertTrue((self.mask.mask == self.test_mask).all())

    def test_getitem(self):
        self.mask['new mask'] = self.test_mask
        self.assertTrue((self.mask['new mask'].mask == self.test_mask).all())
