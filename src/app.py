from flask import Flask
from flask_restful import Api

from src.resources.entity import Entity
from src.resources.event import EventEndpoint
from src.resources.blob import Blob
from src.resources.libamp import Libamp
from src.resources.guide import GuideEndpoint
from src.resources.wge import WGEEndpoint
from src.resources.targeton_oligos import TargetonOligoEndpoint

app = Flask(__name__)
api = Api(app)

@app.route("/")
def hello():
    return "SGE Restful Layer"

api.add_resource(Entity, '/entity/<string:id>')
api.add_resource(EventEndpoint, '/event', methods=["POST"])
api.add_resource(Blob, '/blob/<string:id>')
api.add_resource(Libamp, '/libamp', methods=["POST"])
api.add_resource(GuideEndpoint, '/guide', methods=["POST"])
api.add_resource(WGEEndpoint, '/wge', methods=["POST", "GET"])
api.add_resource(TargetonOligoEndpoint, '/targeton-oligo', methods=["POST"])

if __name__ == "__main__":
    app.run()
