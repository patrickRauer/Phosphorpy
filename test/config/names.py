import unittest

from Phosphorpy.config import names


class TestNames(unittest.TestCase):

    def test_get_available_surveys(self):
        surveys = names.get_available_surveys()
        self.assertIsNotNone(surveys)

        self.assertGreater(len(surveys), 0)

    def test_get_survey_magnitude_names(self):

        for s in names.get_available_surveys():
            cols = names.get_survey_magnitude_names(s)
            self.assertIsNotNone(cols)

        with self.assertRaises(AttributeError):
            names.get_survey_magnitude_names('TEST')
