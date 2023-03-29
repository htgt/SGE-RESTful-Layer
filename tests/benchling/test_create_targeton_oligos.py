import unittest

from src.biology.targeton_oligos import TargetonOligo
from src.benchling.targeton_oligos import (
    as_benchling_entity,
)


class TestTargetonOligo(unittest.TestCase):
    def setUp(self):
        self.data = {
            'entity_name': 'chr17_43104794_43105038_minus_sgRNA_ex4',
            'chromosome': 'chr17',
            'strand': '-',
            'start': '43104794',
            'end': '43105038',
            'r2_start': '43104868',
            'r2_end': '43104956',
            'ext_vector': '25,25',
            'action_vector': '(1del,2del1,snv),(1del,snvre,inframe,stop,ala),(1del,2del0,snv)',
            'sgrna_vector': 'sgRNA_ex4'
        }
        self.targeton_oligo = TargetonOligo(self.data)


    def test_as_benchling_entity_success(self):
        expected_response = {
            'Action Vector': {
                'value': '(1del,2del1,snv),(1del,snvre,inframe,stop,ala),(1del,2del0,snv)'
            },
            'Ext Vector': {
                'value': '25,25'
            },
            'Ref. Start Position': {
                'value': 43104794
            },
            'Ref. End Position': {
                'value': 43105038
            },
            'R2 Start Position': {
                'value': 43104868
            },
            'R2 End Position': {
                'value': 43104956
            },
            'Reference Chromosome': {
                'value': 'sfso_yMb5PNqX'
            },
            'Reference Strand': {
                'value': 'sfso_qKNl7o1M'
            },
            'sgRNA Vector': {
                'value': 'sgRNA_ex4'
            }
        }

        actual_response = as_benchling_entity(self.targeton_oligo)

        self.assertDictEqual(actual_response, expected_response)


if __name__ == '__main__':
    unittest.main()
