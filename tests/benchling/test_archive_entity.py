import unittest
from unittest.mock import patch
from src.benchling.archive_entity import (
    archive_entity, 
    archive_oligo, 
    send_archive_request, 
    prepare_archive_json, 
    generate_entity_ids_str
)
from tests import benchling_ids
from tests.mock_connection import MockRequests, MockResponse


class TestArchiveOligo(unittest.TestCase):
    def setUp(self):
        self.example_type = "dnaOligo"
        self.example_id = benchling_ids['schemas']['grna_oligo_schema_id']
        self.example_entity_ids_str = self.example_type + "Ids"
        self.example_json = {self.example_entity_ids_str: [self.example_id], "reason": "Made in error"}
        self.example_response_json = {
            'dnaOligoIds': [self.example_id], 
            'reason': 'Made in error'
        }
        self.example_status_code = 200
        self.example_status_code_fail = 404
        self.example_url = 'http://test.com'
        self.example_token = 'token'
        self.example_response = MockResponse(self.example_response_json, self.example_status_code)
        self.example_response_fail = MockResponse(self.example_response_json, self.example_status_code_fail)
        self.example_response_text = self.example_response.text
    
    # @patch('src.benchling.archive_entity.archive_entity')
    @patch('src.benchling.archive_entity.archive_entity')
    @patch('src.benchling.connection.connection_class.BenchlingConnection')
    def test_archive_oligo(self, mocked_connection, mocked_archive_entity):
        # Arrange
        example_id = self.example_id
        expected_result = self.example_response_text
        mocked_archive_entity.return_value = expected_result
        # Act
        actual_result = archive_oligo(example_id, mocked_connection)
        # Assert
        self.assertEqual(actual_result, expected_result)
        assert mocked_archive_entity.called
        
    @patch('src.benchling.archive_entity.send_archive_request')
    def test_archive_entity(self, mocked_send_archive_request):
        # Arrange
        example_id = self.example_id
        example_type = self.example_type
        expected_result = self.example_response_text
        mocked_send_archive_request.return_value = expected_result
        example_url = self.example_url
        example_token = self.example_token
        # Act
        actual_result = archive_entity(example_id, example_type, example_url, example_token)
        # Assert
        self.assertEqual(actual_result, expected_result)
        assert mocked_send_archive_request.called
    
    @patch('requests.post', side_effect = MockRequests.post)
    @patch('src.benchling.connection.connection_class.BenchlingConnection')
    def test_send_archive_request(self, mocked_connection, mocked_request_post):
        # Arrange
        example_url = self.example_url
        example_json = self.example_json
        example_token = self.example_token
        expected_result = self.example_response_json
        expected_status_code = self.example_status_code
        # Act
        actual_result = send_archive_request(example_url, example_json, example_token)
        # Assert
        self.assertEqual(actual_result.json(), expected_result)
        self.assertEqual(actual_result.status_code, expected_status_code)
        assert mocked_request_post.called
    
    def test_prepare_archive_json(self):
        # Arrange
        example_type = self.example_type
        example_id = self.example_id
        expected_json = self.example_json
        # Act
        actual_json = prepare_archive_json(example_type, example_id)
        # Assert
        self.assertDictEqual(actual_json, expected_json)
            
    def test_generate_entity_ids_str(self):
        # Arrange
        example_type = self.example_type
        expected_str = self.example_entity_ids_str
        # Act
        actual_str = generate_entity_ids_str(example_type)
        # Assert
        self.assertEqual(actual_str, expected_str)


if __name__ == '__main__':
    unittest.main()
