from flask_restful import Resource

class Entity(Resource):
    def get(self, id):
        return {'entity': id}