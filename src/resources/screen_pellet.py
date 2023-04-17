from flask import request
from flask_restful import Resource

from src.benchling.patch_screen_pellet import patch_screen_pellet

import requests
import json

DATA = {
    "sample_supplier_id": "bfi_nVtf98eI",
    "irods_data_relative_path": "www",
    "run_status": 20,
    "id_run": 1,
}


class ScreenPelletEndpoint(Resource):
    def post(self):
        data = request.json

        response = patch_screen_pellet(data, data["sample_supplier_id"])

        return response

    def get(self):
        patch_screen_pellet(DATA, DATA["sample_supplier_id"])

        return DATA
