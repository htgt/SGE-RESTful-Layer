from flask_restful import Resource
import json

from src.benchling.get_blob import get_blob_url
from src.benchling.guideRNA_from_csv import GrnasImportFromCSV


class Blob(Resource):
    def get(self, id):
        csv_url = get_blob_url(id)

        result = GrnasImportFromCSV().get_grnas(csv_url)

        return result, 201
