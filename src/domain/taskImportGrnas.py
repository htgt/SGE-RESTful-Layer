from src.benchling.get_blob import get_blob_url


class AvailableTaskStatuses:
    def __init__(self, pending_id, failed_id, completed_id):
        self.pending = pending_id
        self.failed = failed_id
        self.completed = completed_id

class TaskImport:
    def __init__(self, data, get_url_method=get_blob_url):
        self.id = data["id"]
        self.status_id = data["status_id"]
        self.file_id = data["file_id"]

        if self.file_id != None:
            try:
                self.file_url = get_url_method(self.file_id)
            except Exception as err:
                raise Exception("Could not get input file url")

    def execute(self):
        try:
            importer = GrnasImportFromCSV()
            result = importer.import_grnas(self.file_url)
        except Exception as err:
            return "Could not import guide RNAs", 500

        return result, 200

    def update_status(self, status):
        pass


