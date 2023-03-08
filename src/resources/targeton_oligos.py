from flask import request
from flask_restful import Resource

from src.benchling.targeton_oligos import post_targeton_oligos

import requests
import json

class TargetonOligoEndpoint(Resource):
    def post(self):
        data = request.json

        post_targeton_oligos(data)

        return data
