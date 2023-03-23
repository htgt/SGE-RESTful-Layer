from flask import request
from flask_restful import Resource
from src.core.guide import handle_guide_event

import json


class GuideEndpoint(Resource):
    def get(self, id):
        return id, 201

    def post(self):
        data = request.json

        return handle_guide_event(data)
