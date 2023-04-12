from flask import request
from flask_restful import Resource
from src.core.guide import handle_guide_event
from src.benchling.connection.benchling_connection import benchling_connection, benchling_schema_ids

import json


class GuideEndpoint(Resource):
    def get(self, id):
        return id, 201

    def post(self):
        data = request.json

        return handle_guide_event(data, benchling_connection, benchling_schema_ids)
