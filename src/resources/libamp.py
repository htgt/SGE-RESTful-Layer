from flask import request
from flask_restful import Resource

class Libamp(Resource):
    def post(self):
        data = request.json

        print(data)

        return data, 201