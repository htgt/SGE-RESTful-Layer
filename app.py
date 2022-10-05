from flask import Flask
from flask_restful import Api

from resources.entity import Entity

app = Flask(__name__)
api = Api(app)

api.add_resource(Entity, '/entity/<string:id>')

app.run()