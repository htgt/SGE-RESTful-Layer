import unittest

from src.domain.taskImportGrnas import TaskImport


class TestRaskImportGrnas(unittest.TestCase):
    def _get_blob_url_fake(self, str):
        return "www.example.com"

    def test_create_task_import(self):
        task_data = {
            "id": "1168686327",
            "status_id": "wfts_DrqRAcOl",
            "file_id"     : "A1BG"
        }

        test_task = TaskImport(task_data, self._get_blob_url_fake)

        self.assertEqual(getattr(test_task, "id"), "1168686327")
        self.assertEqual(getattr(test_task, "status_id"),
                         'wfts_DrqRAcOl')
        self.assertEqual(getattr(test_task, "file_id"), "A1BG")
        self.assertEqual(getattr(test_task, "file_url"), "www.example.com")
