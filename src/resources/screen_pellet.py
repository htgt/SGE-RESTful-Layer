from flask import request
from flask_restful import Resource

from src.benchling.patch_screen_pellet import patch_screen_pellet
from src.benchling.connection.benchling_connection import benchling_connection

import requests
import json


class ScreenPelletEndpoint(Resource):
    def post(self):
        data = request.json

        response = patch_screen_pellet(data, data["sample_supplier_id"], benchling_connection)

        return response
