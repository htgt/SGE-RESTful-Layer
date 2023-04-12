from flask_restful import Resource
import json

from src.benchling.get_blob import get_blob_url
from src.benchling.guideRNA_from_csv import GrnasImportFromCSV
from src.benchling.connection.benchling_connection import benchling_connection


class Blob(Resource):
    def get(self, id):
        api_path = benchling_connection.blobs_url
        token = benchling_connection.token
        csv_url = get_blob_url(id, api_path, token)

        result = GrnasImportFromCSV().get_grnas(csv_url)

        return result, 201
