from flask import request
from flask_restful import Resource
from src.core.guide import handle_guide_event
from src.benchling.connection.benchling_connection import benchling_connection

import json


class GuideEndpoint(Resource):
    def get(self, id):
        return id, 201

    def post(self):
        data = request.json
        url = benchling_connection.oligos_url
        return handle_guide_event(data, url)
