from flask import request
from flask_restful import Resource
from src.core.guide import patch_grna_event, post_grna_oligos_event

import json

BENCHLING_GUIDE_RNA_SCHEMA_ID = "ts_vGZYroiQ"
BENCHLING_ENTITY_REGISTERED_EVENT = "v2.entity.registered"


class GuideEndpoint(Resource):
    def get(self, id):

        return id, 201

    def post(self):
        data = request.json
        response = {}
        
        try:
            response['grna'] = patch_grna_event(data)
        except Exception as err:
            return json.dumps(err), 500
        try:
            response['oligos'] = post_grna_oligos_event(data)
            return response, 200
        except Exception as err:
            return json.dumps(err), 500
