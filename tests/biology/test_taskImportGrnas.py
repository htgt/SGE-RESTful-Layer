import unittest
from unittest.mock import patch
from src.biology.taskImportGrnas import TaskImport


class TestRaskImportGrnas(unittest.TestCase):
    def _get_blob_url_fake(self, *args):
        return "www.example.com"

    def test_create_task_import(self):
        task_data = {
            "id": "1168686327",
            "status_id": "wfts_DrqRAcOl",
            "file_id"     : "A1BG"
        }
        api_path = 'www.fake_api.test'
        token = 'test'
        test_task = TaskImport(task_data, api_path, token, self._get_blob_url_fake)

        self.assertEqual(getattr(test_task, "id"), "1168686327")
        self.assertEqual(getattr(test_task, "status_id"),
                         'wfts_DrqRAcOl')
        self.assertEqual(getattr(test_task, "file_id"), "A1BG")
        self.assertEqual(getattr(test_task, "file_url"), "www.example.com")
