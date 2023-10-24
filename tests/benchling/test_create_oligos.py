import unittest
from unittest.mock import patch
from src.utils.exceptions import OligoDirectionInvalid
from src.biology.guideRNA import Oligo, GuideRNAOligos
from src.benchling.create_oligos import prepare_oligo_json, setup_oligo_class, setup_oligo_pair_class, export_oligos_to_benchling, BenchlingOligo, BenchlingOligosPair 
from Bio.Seq import Seq
from copy import deepcopy
from tests import benchling_ids
from tests.mock_connection import Mock_request_to_benchling


class TestCreateOligo(unittest.TestCase):
    def setUp(self):
        self.example_oligo_schema_id = benchling_ids['schemas']['grna_oligo_schema_id']
        self.example_grna_schema_id = benchling_ids['schemas']['grna_schema_id']
        self.example_sense_id = benchling_ids['dropdowns']['sense']
        self.example_antisense_id = benchling_ids['dropdowns']['antisense']
        self.example_targeton_id = benchling_ids['schemas']['targeton_oligo_schema_id']
        self.example_folder_id = 'lib_MaKCkHDE'
        
        self.example_seq = 'AGCTGACGGGTGACACCCC'
        self.example_seq_forward = 'CACCGGCTGACGGGTGACACCCC'
        self.example_seq_reverse = 'AAACGGGGTGTCACCCGTCAGCC'
        self.example_guide_seq = 'TGCTGACGGGTGACACCCA'
        
        self.example_forward_oligo_json_dict = {
            'bases': self.example_seq_forward,
            'fields': {
                'Targeton': {'value': self.example_targeton_id},
                'Strand': {'value': self.example_sense_id},
                'Guide RNA': {'value': self.example_grna_schema_id}
            },
            'folderId': self.example_folder_id,
            'name': 'Guide RNA Oligo',
            'schemaId': self.example_oligo_schema_id,
            'isCircular': False,
        }
        self.example_reverse_oligo_json_dict = {
            'bases': self.example_seq_reverse,
            'fields': {
                'Targeton': {'value': self.example_targeton_id},
                'Strand': {'value': self.example_antisense_id},
                'Guide RNA': {'value': self.example_grna_schema_id}
            },
            'folderId': self.example_folder_id,
            'name': 'Guide RNA Oligo',
            'schemaId': self.example_oligo_schema_id,
            'isCircular': False,
        }
        self.example_oligos_pair = GuideRNAOligos(self.example_seq)
        self.example_benchling_oligos_pair = BenchlingOligosPair(
            forward=BenchlingOligo(
                sequence=Seq(self.example_seq_forward),
                targeton=self.example_targeton_id,
                folder_id=self.example_folder_id,
                schema_id=self.example_oligo_schema_id,
                name='Guide RNA Oligo',
                strand=self.example_sense_id,
                grna=self.example_grna_schema_id
            ),
            reverse=BenchlingOligo(
                sequence=Seq(self.example_seq_reverse),
                targeton=self.example_targeton_id,
                folder_id=self.example_folder_id,
                schema_id=self.example_oligo_schema_id,
                name='Guide RNA Oligo',
                strand=self.example_antisense_id,
                grna=self.example_grna_schema_id
            )
        )
        self.example_guide_data = {
            'id': self.example_grna_schema_id,
            'targeton': self.example_targeton_id,
            'folder_id': self.example_folder_id,
            'schemaid': self.example_oligo_schema_id,
            'name': 'Guide RNA Oligo',
            'seq': self.example_guide_seq
        }
        self.example_oligo = Oligo(sequence=Seq(self.example_seq_forward))
        self.example_oligo_forward_dict = {
                'sequence': Seq(self.example_seq_forward),
                'targeton': self.example_targeton_id,
                'folder_id': self.example_folder_id,
                'schema_id': self.example_oligo_schema_id,
                'name': 'Guide RNA Oligo',
                'strand': self.example_sense_id,
                'grna': self.example_grna_schema_id,
        }
        self.example_oligo_reverse_dict = {
                'sequence': Seq(self.example_seq_reverse),
                'targeton': self.example_targeton_id,
                'folder_id': self.example_folder_id,
                'schema_id': self.example_oligo_schema_id,
                'name': 'Guide RNA Oligo',
                'strand': self.example_antisense_id,
                'grna': self.example_grna_schema_id,
        }
        self.example_setup_oligos_list_dicts = [
            self.example_oligo_forward_dict,
            self.example_oligo_reverse_dict
        ]

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
        expected_forward = self.example_oligo_forward_dict
        expected_reverse = self.example_oligo_reverse_dict
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
        actual_forward = forward._asdict()
        actual_reverse = reverse._asdict()
        # Assert
        self.assertDictEqual(actual_forward, expected_forward)
        self.assertDictEqual(actual_reverse, expected_reverse)

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
    
    @patch.dict('src.benchling.benchling_schema_ids.ids', benchling_ids)
    def test_setup_oligo_pair_class(self):
        # Arrange
        example_oligos = self.example_oligos_pair
        example_guide_data = self.example_guide_data
        expected_result = self.example_setup_oligos_list_dicts
        # Act
        actual_result = setup_oligo_pair_class(example_oligos, example_guide_data)
        # Assert 
        self.assertCountEqual(actual_result.to_list_dicts(), expected_result)
    
    @patch('src.benchling.utils.request_to_benchling.request_to_benchling', side_effect = Mock_request_to_benchling)
    def test_export_oligos_to_benchling(self, mocked_request_benchling):
        # Arrange
        example_oligos = self.example_benchling_oligos_pair
        expected_response = (
            {'bases': self.example_seq_forward, 'fields': {'Targeton': {'value': self.example_targeton_id}, 'Strand': {'value': self.example_sense_id}, 'Guide RNA': {'value': self.example_grna_schema_id}}, 'isCircular': False, 'folderId': self.example_folder_id, 'name': 'Guide RNA Oligo', 'schemaId': self.example_oligo_schema_id},
            {'bases': self.example_seq_reverse, 'fields': {'Targeton': {'value': self.example_targeton_id}, 'Strand': {'value': self.example_antisense_id}, 'Guide RNA': {'value': self.example_grna_schema_id}}, 'isCircular': False, 'folderId': self.example_folder_id, 'name': 'Guide RNA Oligo', 'schemaId': self.example_oligo_schema_id}
            )
        # Act
        actual_response = export_oligos_to_benchling(example_oligos)
        # Assert
        self.assertCountEqual(actual_response, expected_response)
        mocked_request_benchling.assert_called()


if __name__ == '__main__':
    unittest.main()
