import unittest
from unittest.mock import patch
from src.utils.exceptions import OligoDirectionInvalid
from src.biology.guideRNA import Oligo, GuideRNAOligos
from src.benchling.create_oligos import prepare_oligo_json, setup_oligo_class, BenchlingOligo, BenchlingOligosPair
import json
from Bio.Seq import Seq
from copy import deepcopy
from tests import benchling_ids


class TestCreateOligo(unittest.TestCase):
    def setUp(self):
        self.example_oligo_schema_id = benchling_ids['schemas']['grna_oligo_schema_id']
        self.example_grna_schema_id = benchling_ids['schemas']['grna_schema_id']
        self.example_sense_id = benchling_ids['dropdowns']['sense']
        self.example_antisense_id = benchling_ids['dropdowns']['antisense']
        
        self.example_forward_oligo_json_dict = {
            'bases': 'CACCGGCTGACGGGTGACACCCC',
            'fields': {
                'Targeton': {'value': 'seq_8VA7PA1S'},
                'Strand': {'value': self.example_sense_id},
                'Guide RNA': {'value': self.example_grna_schema_id}
            },
            'folderId': 'lib_MaKCkHDE',
            'name': 'Guide RNA Oligo',
            'schemaId': self.example_oligo_schema_id,
            'isCircular': False,
        }
        self.example_reverse_oligo_json_dict = {
            'bases': 'AAACGGGGTGTCACCCGTCAGCC',
            'fields': {
                'Targeton': {'value': 'seq_8VA7PA1S'},
                'Strand': {'value': self.example_antisense_id},
                'Guide RNA': {'value': self.example_grna_schema_id}
            },
            'folderId': 'lib_MaKCkHDE',
            'name': 'Guide RNA Oligo',
            'schemaId': self.example_oligo_schema_id,
            'isCircular': False,
        }
        
        self.example_seq = 'AGCTGACGGGTGACACCCC'
        self.example_oligos_pair = GuideRNAOligos(self.example_seq)
        self.example_benchling_oligos_pair = BenchlingOligosPair(
            forward=BenchlingOligo(
                sequence=Seq('CACCGGCTGACGGGTGACACCCC'),
                targeton='seq_8VA7PA1S',
                folder_id='lib_MaKCkHDE',
                schema_id=self.example_oligo_schema_id,
                name='Guide RNA Oligo',
                strand=self.example_sense_id,
                grna=self.example_grna_schema_id
            ),
            reverse=BenchlingOligo(
                sequence=Seq('AAACGGGGTGTCACCCGTCAGCC'),
                targeton='seq_8VA7PA1S',
                folder_id='lib_MaKCkHDE',
                schema_id=self.example_oligo_schema_id,
                name='Guide RNA Oligo',
                strand=self.example_antisense_id,
                grna=self.example_grna_schema_id
            )
        )
        self.example_guide_data = {
            'id': self.example_grna_schema_id,
            'targeton': 'seq_8VA7PA1S',
            'folder_id': 'lib_MaKCkHDE',
            'schemaid': self.example_oligo_schema_id,
            'name': 'Guide RNA Oligo',
            'seq': 'TGCTGACGGGTGACACCCA'
        }
        self.example_oligo = Oligo(sequence=Seq('CACCGGCTGACGGGTGACACCCC'))
        self.example_setup_oligos_list_dicts = [
            {
                'sequence': Seq('CACCGGCTGACGGGTGACACCCC'),
                'targeton': 'seq_8VA7PA1S',
                'folder_id': 'lib_MaKCkHDE',
                'schema_id': self.example_oligo_schema_id,
                'name': 'Guide RNA Oligo',
                'strand': self.example_sense_id,
                'grna': self.example_grna_schema_id,
            },
            {
                'sequence': Seq('AAACGGGGTGTCACCCGTCAGCC'),
                'targeton': 'seq_8VA7PA1S',
                'folder_id': 'lib_MaKCkHDE',
                'schema_id': self.example_oligo_schema_id,
                'name': 'Guide RNA Oligo',
                'strand': self.example_antisense_id,
                'grna': self.example_grna_schema_id,
            }
        ]

    # This method will be used by the mock to export_to_benchling
    def mocked_requests_get(*args, **kwargs):
        class MockResponse:
            def __init__(self, json_data, status_code):
                self.json_data = json_data
                self.status_code = status_code

            def json(self):
                return self.json_data

        if args[0] == 'http://someurl.com/test.json':
            return MockResponse({"key1": "value1"}, 200)
        elif args[0] == 'http://someotherurl.com/anothertest.json':
            return MockResponse({"key2": "value2"}, 200)

        return MockResponse(None, 404)

    def test_prepare_oligo_json(self):
        # Arrange
        oligos = deepcopy(self.example_benchling_oligos_pair)
        # Act
        test_forward_oligo_json_dict = prepare_oligo_json(oligos.forward)
        test_reverse_oligo_json_dict = prepare_oligo_json(oligos.reverse)
        # Assert
        example_forward_oligo_json_dict = self.example_forward_oligo_json_dict
        example_reverse_oligo_json_dict = self.example_reverse_oligo_json_dict

        self.assertDictEqual(test_forward_oligo_json_dict, example_forward_oligo_json_dict)
        self.assertDictEqual(test_reverse_oligo_json_dict, example_reverse_oligo_json_dict)

    def test_setup_oligo_class(self):
        # Arrange
        oligos = deepcopy(self.example_oligos_pair)
        guide_data = self.example_guide_data
        # Act
        # Foward
        forward = setup_oligo_class(
            oligos.forward,
            guide_data,
            benchling_ids,
            'forward',
        )
        # Reverse
        reverse = setup_oligo_class(
            oligos.reverse,
            guide_data,
            benchling_ids,
            'reverse',
        )
        benchling_oligos = BenchlingOligosPair(forward, reverse)
        # Assert
        test_oligos_list_dicts = benchling_oligos.to_list_dicts()
        example_oligos_list_dicts = self.example_setup_oligos_list_dicts
        
        self.assertCountEqual(test_oligos_list_dicts, example_oligos_list_dicts)

    def test_wrong_direction_setup_oligo_class(self):
        # Arrange
        oligos = deepcopy(self.example_oligos_pair)
        guide_data = self.example_guide_data
        # Act, Assert
        with self.assertRaises(OligoDirectionInvalid):
            oligos.forward = setup_oligo_class(
                oligos.forward,
                guide_data,
                benchling_ids,
                'wrong_direction',
            )


if __name__ == '__main__':
    unittest.main()
