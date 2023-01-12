from flask import request
from flask_restful import Resource

entities = []


class Entity(Resource):
    def get(self, id):
        entity = None

        for item in entities:
            if item['id'] == id:
                entity = item

        return {'entity': entity}, 200 if entity else 404

    def post(self, id):
        data = request.get_json()
        entity = {'id': id, 'name': data['name']}
        entities.append(entity)

        return entity, 201
