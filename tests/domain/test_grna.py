import unittest
import pdb 

from Bio.Seq import Seq
from src.domain.guideRNA import GuideRNA, GuideRNAOligo, create_set_of_gRNAs


class TestGuideRNA(unittest.TestCase):
    def test_create_guide_RNA(self):
        input_data = {
            'wge_id': '1168686327',
            'seq': 'GACTTCCAGCTACGGCGCG',
            'targeton': 'TGTN001',
            'strand': 'sfso_qKNl7o1M',
            'wge_link': 'www.test.com',
            'off_targets': '{0: 1, 1: 0, 2: 1, 3: 15, 4: 204}',
            'species': 'sfso_gWKuC1ge',
        }

        test_gRNA = GuideRNA(input_data)

        self.assertEqual(getattr(test_gRNA, "wge_id"), "1168686327")
        self.assertEqual(getattr(test_gRNA, "sequence"),
                         Seq("GACTTCCAGCTACGGCGCG"))
        self.assertEqual(getattr(test_gRNA, "off_targets"), "{0: 1, 1: 0, 2: 1, 3: 15, 4: 204}")

    def test_grna_as_benchling_req_body(self):
        input_data = {
            'wge_id': '1168686327',
            'seq': 'GACTTCCAGCTACGGCGCG',
            'targeton': 'TGTN001',
            'strand': 'sfso_qKNl7o1M',
            'wge_link': 'www.test.com',
            'off_targets': '{0: 1, 1: 0, 2: 1, 3: 15, 4: 204}',
            'species': 'sfso_gWKuC1ge',
        }
        
        event = {
            'folder_id' : 'folder',
            'name' : 'test',
            'schema_id' : 'schema_id',
        }

        expected_out = {
            'bases': 'GACTTCCAGCTACGGCGCG',
            'fields': {
                'WGE ID': {
                    'value': '1168686327'
                },
                'Guide Sequence': {
                    'value': 'GACTTCCAGCTACGGCGCG'
                },
                'Targeton': {
                    'value': 'TGTN001'
                },
                'Strand': {
                    'value': 'sfso_qKNl7o1M'
                },
                'WGE Hyperlink': {
                    'value': 'www.test.com'
                },
                'Off Target Summary Data': {
                    'value': '{0: 1, 1: 0, 2: 1, 3: 15, 4: 204}'
                },
                'Species': {
                    'value': 'sfso_gWKuC1ge'
                }
            },
            'folderId': 'folder',
            'name': 'test',
            'schemaId': 'schema_id'
        }

        test_gRNA = GuideRNA(input_data)
        
        self.assertEqual(test_gRNA.as_benchling_req_body(event), expected_out)


class TestGuideRNAOligo(unittest.TestCase):
    def test_forward_sequence(self):
        input_sequence = 'AATATGGTGGCCCTCCATT'
        first_base = 'G'

        transformed = GuideRNAOligo(input_sequence).forward_sequence()

        self.assertEqual(transformed[0], first_base)
        self.assertEqual(transformed, 'GATATGGTGGCCCTCCATT')

    def test_reverse_sequence(self):
        input_sequence = 'AATATGGTGGCCCTCCATT'
        last_base = 'C'

        transformed = GuideRNAOligo(input_sequence).reverse_sequence()

        self.assertEqual(transformed[-1], last_base)
        self.assertEqual(transformed, 'AATGGAGGGCCACCATATC')

    def test_create_oligos(self):
        input_sequence = 'AATATGGTGGCCCTCCATT'

        oligos = GuideRNAOligo(input_sequence).create_oligos()

        self.assertEqual(oligos.forward.sequence, Seq('CACCGATATGGTGGCCCTCCATT'))
        self.assertEqual(oligos.reverse.sequence, Seq('AAACAATGGAGGGCCACCATATC'))


if __name__ == '__main__':
    unittest.main()
