from flask import request
from flask_restful import Resource

from src.benchling.post_libamp_primers import post_libamp_primers

class Libamp(Resource):
    def post(self):
        data = request.json

        #print(data)

        result = post_libamp_primers(data)

        return result, 201