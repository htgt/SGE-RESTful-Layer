from flask_restful import Resource
from src.benchling.get_blob import get_blob_url


class Blob(Resource):
    def get(self, id):
        url = get_blob_url(id)

        return url, 201
