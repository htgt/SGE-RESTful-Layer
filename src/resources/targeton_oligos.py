from flask import request
from flask_restful import Resource

from src.benchling.connection.benchling_connection import benchling_connection
from src.benchling.targeton_oligos import post_targeton_oligos


class TargetonOligoEndpoint(Resource):
    def post(self):
        data = request.json

        post_targeton_oligos(data, benchling_connection.custom_entity_url)

        return data
