import unittest

from Phosphorpy.config import survey_data


class TestConfig(unittest.TestCase):

    def test_create_dict(self):
        self.assertEqual(survey_data.create_dict('test 1'),
                         {'test': '1'})
        self.assertEqual(survey_data.create_dict('test   1'),
                         {'test': '1'})
        self.assertEqual(survey_data.create_dict('test 1 2 3'),
                         {'test': ['1', '2', '3']})

    def test_read_survey_data(self):
        self.assertIsNotNone(survey_data.read_survey_data())
