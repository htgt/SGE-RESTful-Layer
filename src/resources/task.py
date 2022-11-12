from flask import request
from flask_restful import Resource

class TaskEndpoint(Resource):
    def get(self, id):
        return id, 201

    def post(self):
        data = request.json

        return data, 201