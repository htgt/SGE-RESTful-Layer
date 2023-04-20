import unittest
from unittest.mock import patch

from src.benchling.patch_screen_pellet import (
    patch_screen_pellet,
    as_benchling_req_body,
    convert_status_id_to_status
)


class TestPatchScreenPellet(unittest.TestCase):
    def setUp(self):
        self.data = {
            "sample_supplier_id": "jkl123",
            "irods_data_relative_path": "path",
            "run_status": 20,
            "id_run": 1,
        }

    @patch('src.benchling.patch_screen_pellet.as_benchling_req_body')
    @patch('src.benchling.patch_screen_pellet.request_to_benchling_json_response')
    def test_patch_screen_pellet(self, mock_export, mock_body):
        # arrange
        mock_export.return_value = 'test response'
        mock_body.return_value = {'test': 'body'}
        url = 'url'
        expected = 'test response'

        # act
        actual = patch_screen_pellet(self.data, 'jkl123', url)

        # assert
        self.assertEqual(actual, expected)
        mock_export.assert_called_with({'test': 'body'}, 'url/jkl123', 'patch')

    @patch('src.benchling.patch_screen_pellet.benchling_schema_ids')
    def test_as_benchling_req_body(self, mock_schema_ids):
        # arrange
        mock_schema_ids.ids = {
            'schemas': {
                'screen_pellet_schema_id': 'abc123',
            },
            'dropdowns': {
                'sequencing_status': {
                    'submitted': 'def456',
                    'finished': 'ghi789',
                },
            },
        }
        expected = {
            'fields': {
                'iRODs Path': {
                    'value': 'path',
                },
                'Sequencing Status': {
                    'value': 'ghi789',
                },
                'Run ID': {
                    'value': 1,
                },
            },
            'schemaId': 'abc123',
        }

        # act
        actual = as_benchling_req_body(self.data)

        # assert
        self.assertEqual(actual, expected)

    def test_convert_status_id_to_status_submitted(self):
        # arrange
        expected = 'submitted'

        # act
        actual = convert_status_id_to_status(5)

        # assert
        self.assertEqual(actual, expected)

    def test_convert_status_id_to_status_finished(self):
        # arrange
        expected = 'finished'

        # act
        actual = convert_status_id_to_status(20)

        # assert
        self.assertEqual(actual, expected)
