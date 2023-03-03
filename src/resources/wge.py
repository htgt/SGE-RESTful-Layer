from flask import request
from flask_restful import Resource

from src.wge.wge import patch_wge_data_to_service

import requests
import json

class WGEEndpoint(Resource):
    def __transform_event_input_data(self, data):
        data_entity = data['detail']['entity']
        
        wge_grna_data = {}
        wge_grna_data['folder_id'] = data_entity['folderId']
        wge_grna_data['entity_id'] = data_entity['id']
        wge_grna_data['wge_id'] = data_entity['fields']['WGE ID']['value']
        wge_grna_data['targeton_id'] = data_entity['fields']['Targeton']['value']
        wge_grna_data['schema_id'] = data_entity['schema']['id']
        wge_grna_data['name'] = data_entity['schema']['name']
        
        print(wge_grna_data)

        return wge_grna_data

    def get(self):
        wge_id = request.args.get('id')
        wge_json = query_wge_by_id(wge_id)

        return wge_json
    
    def post(self):
        data = request.json
        event_data = self.__transform_event_input_data(data)
        response = patch_wge_data_to_service(self, event_data)

        return response
