from src.benchling.get_blob import get_blob_url
from src.benchling.guideRNA_from_csv import GrnasImportFromCSV

class TaskImport:
    def __init__(self, data):
        self.id = data["id"]
        self.file_id = data["file_id"]
        self.status_id = data["status_id"]

    def execute(self):
        try:
            file_url = get_blob_url(self.file_id)

            importer = GrnasImportFromCSV()
            result = importer.import_grnas(file_url)
        except Exception as err:
            return "Could not import guide RNAs", 500

        return result, 200

    def update_status(self, status):
        pass