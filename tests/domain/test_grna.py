import unittest
import pdb 

from Bio.Seq import Seq
from src.biology.guideRNA import GuideRNA, GuideRNAOligos
from src.benchling.patch_guide_rna import as_benchling_req_body
import json


class TestGuideRNA(unittest.TestCase):
    def setUp(self):
        with open('benchling_schema_ids.json', 'r') as f:
            self.benchling_ids = json.load(f)
        self.example_seq = 'GACTTCCAGCTACGGCGCG'
        self.example_wge_id = '1168686327'
        self.example_wge_link = 'www.test.com'
        self.example_off_targets = '{0: 1, 1: 0, 2: 1, 3: 15, 4: 204}',
        self.example_species_benchling_id = self.benchling_ids['dropdowns']['species']['homo_sapiens']
        self.example_species = 'Grch37'
        self.example_strand = 'sfso_qKNl7o1M'
        self.example_targeton = 'TGTN001'
        self.example_input_data = {
            'wge_id': self.example_wge_id,
            'seq': self.example_seq,
            'targeton': self.example_targeton,
            'strand': self.example_strand,
            'wge_link': self.example_wge_link,
            'off_targets': self.example_off_targets,
            'species': self.example_species,
        }
    
    def test_create_guide_RNA(self):
        input_data = self.example_input_data

        test_gRNA = GuideRNA(input_data)

        self.assertEqual(getattr(test_gRNA, "wge_id"), "1168686327")
        self.assertEqual(getattr(test_gRNA, "sequence"),
                         Seq("GACTTCCAGCTACGGCGCG"))
        self.assertEqual(getattr(test_gRNA, "off_targets"), "{0: 1, 1: 0, 2: 1, 3: 15, 4: 204}")

    def test_grna_as_benchling_req_body(self):
        # Arrange
        input_data = self.example_input_data
        
        event = {
            'folder_id' : 'folder',
            'name' : 'test',
            'schema_id' : 'schema_id',
            'targeton_id' : 'targeton_id'
        }

        expected_out = {
            'bases': self.example_seq,
            'fields': {
                'WGE ID': {
                    'value': self.example_wge_id
                },
                'Guide Sequence': {
                    'value': self.example_seq
                },
                'Targeton': {
                    'value': self.example_targeton
                },
                'Strand': {
                    'value': self.example_strand
                },
                'WGE Hyperlink': {
                    'value': self.example_wge_link
                },
                'Off Target Summary Data': {
                    'value': self.example_off_targets
                },
                'Species': {
                    'value': self.example_species
                }
            },
            'folderId': 'folder',
            'name': 'test',
            'schemaId': 'schema_id',
            'Targeton': {'value': 'targeton_id'}
        }
        # Act
        test_gRNA = GuideRNA(input_data)
        # Assert
        self.assertEqual(as_benchling_req_body(test_gRNA, event), expected_out)


class TestGuideRNAOligos(unittest.TestCase):
    def test_forward_sequence(self):
        input_sequence = 'AATATGGTGGCCCTCCATT'
        first_base = 'G'

        transformed = GuideRNAOligos(input_sequence).forward_sequence()

        self.assertEqual(transformed[0], first_base)
        self.assertEqual(transformed, 'GATATGGTGGCCCTCCATT')

    def test_reverse_sequence(self):
        input_sequence = 'AATATGGTGGCCCTCCATT'
        last_base = 'C'

        transformed = GuideRNAOligos(input_sequence).reverse_sequence()

        self.assertEqual(transformed[-1], last_base)
        self.assertEqual(transformed, 'AATGGAGGGCCACCATATC')

    def test_create_oligos(self):
        input_sequence = 'AATATGGTGGCCCTCCATT'

        oligos = GuideRNAOligos(input_sequence)

        self.assertEqual(oligos.forward.sequence, Seq('CACCGATATGGTGGCCCTCCATT'))
        self.assertEqual(oligos.reverse.sequence, Seq('AAACAATGGAGGGCCACCATATC'))


if __name__ == '__main__':
    unittest.main()
