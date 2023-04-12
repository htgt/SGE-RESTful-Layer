import unittest
from src.benchling.patch_guide_rna import as_benchling_req_body
from src.biology.guideRNA import GuideRNA
from src.benchling import BenchlingSchemaIds
import json

class TestPatchGuideRNA(unittest.TestCase):
    def setUp(self):
        self.benchling_ids = BenchlingSchemaIds('tests/fixtures/example_benchling_schema_ids.json')
        self.example_seq = 'GACTTCCAGCTACGGCGCG'
        self.example_wge_id = '1168686327'
        self.example_wge_link = 'www.test.com'
        self.example_off_targets = '{0: 1, 1: 0, 2: 1, 3: 15, 4: 204}',
        self.example_species_benchling_id = self.benchling_ids.ids['dropdowns']['species']['homo_sapiens']
        self.example_species = 'Grch37'
        self.example_strand = 'sfso_qKNl7o1M'
        self.example_targeton = 'TGTN001'
        self.example_input_data = {
            'wge_id': self.example_wge_id,
            'seq': self.example_seq,
            'targeton': self.example_targeton,
        #    'strand': self.example_strand,
            'wge_link': self.example_wge_link,
            'off_targets': self.example_off_targets,
            'species': self.example_species,
        }

    def test_grna_as_benchling_req_body(self):
        # Arrange
        input_data = self.example_input_data

        event = {
            'folder_id': 'folder',
            'name': 'test',
            'schema_id': 'schema_id',
            'targeton_id': 'targeton_id'
        }

        expected_out = {
            'bases': self.example_seq,
            'fields': {
                'WGE ID': {
                    'value': self.example_wge_id
                },
                'Targeton': {
                    'value': self.example_targeton
                },
                'WGE Hyperlink': {
                    'value': self.example_wge_link
                },
                'Off Target Summary Data': {
                    'value': self.example_off_targets
                },
                'Species': {
                    'value': self.example_species_benchling_id
                },
                'Targeton': {
                    'value': 'targeton_id'
                }
            },
            # 'folderId': 'folder',
            'name': 'test',
            'schemaId': 'schema_id',

        }
        # Act
        test_gRNA = GuideRNA(input_data)
        benchling_guide_rna = as_benchling_req_body(test_gRNA, event, self.benchling_ids)

        # Assert
        self.assertEqual(benchling_guide_rna, expected_out)