from src.benchling.get_blob import get_blob_url


class TaskImport:
    def __init__(self, data, api_path, token, get_url_method=get_blob_url):
        self.id = data["id"]
        self.status_id = data["status_id"]
        self.file_id = data["file_id"]

        if self.file_id != None:
            try:
                self.file_url = get_url_method(self.file_id, api_path, token)
            except Exception as err:
                raise Exception(f"Could not get input file url, {err}")

    def get_blob_url(self, get_url_method, *args):
        return get_url_method(*args)

    def execute(self):
        pass

    def update_status(self, status):
        pass
