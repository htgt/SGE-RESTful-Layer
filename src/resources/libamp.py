from flask import request
from flask_restful import Resource
from flask import make_response, jsonify

from src.core.post_libamp_primers import post_libamp_primers

class Libamp(Resource):
    def post(self):
        data = request.json
        result = post_libamp_primers(data)

        return result, 200