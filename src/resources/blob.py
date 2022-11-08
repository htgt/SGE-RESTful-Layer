from flask_restful import Resource
import json

from src.benchling.get_blob import get_blob_url
from src.benchling.guideRNA_from_csv import GrnasImportFromCSV

class Blob(Resource):
    def get(self, id):
        blob = get_blob_url(id)
        url = json.loads(blob)["downloadURL"]

        result = GrnasImportFromCSV().import_grnas(url)

        return 201
