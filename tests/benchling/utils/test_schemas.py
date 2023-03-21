import unittest

from src.benchling.utils.schemas import (
    get_strand_dropdown_id,
    get_chromosome_dropdown_id
)


class TestTargetonOligo(unittest.TestCase):
    def test_get_strand_dropdown_id_plus_success(self):
        input_chr = '+'
        expected_response = 'sfso_DqRsZ1Cg'

        actual_response = get_strand_dropdown_id(input_chr)

        self.assertEqual(actual_response, expected_response)


    def test_get_strand_dropdown_id_minus_success(self):
        input_chr = '-'
        expected_response = 'sfso_qKNl7o1M'

        actual_response = get_strand_dropdown_id(input_chr)

        self.assertEqual(actual_response, expected_response)


    def test_get_chromosome_dropdown_id_chrint_success(self):
        input_chr = 'chr11'
        expected_response = 'sfso_PH8TVTOe'

        actual_response = get_chromosome_dropdown_id(input_chr)

        self.assertEqual(actual_response, expected_response)


    def test_get_chromosome_dropdown_id_chrchr_success(self):
        input_chr = 'chrX'
        expected_response = 'sfso_wVzT3oIa'

        actual_response = get_chromosome_dropdown_id(input_chr)

        self.assertEqual(actual_response, expected_response)
