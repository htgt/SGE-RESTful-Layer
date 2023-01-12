from flask import Flask
from flask_restful import Api

from src.resources.entity import Entity
from src.resources.event import EventEndpoint
from src.resources.blob import Blob
from src.resources.task import TaskEndpoint

app = Flask(__name__)
api = Api(app)

api.add_resource(Entity, '/entity/<string:id>')
api.add_resource(EventEndpoint, '/event', methods=["POST"])
api.add_resource(Blob, '/blob/<string:id>')
api.add_resource(TaskEndpoint, '/task', methods=["POST", "GET"])

if __name__ == "__main__":
    app.run()
