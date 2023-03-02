from flask import request
from flask_restful import Resource
from src.wge.wge import (
    query_wge_by_id,
    prepare_guide_rna_class,
)
from src.rest_calls.send_calls import export_to_service
from src.benchling import benchling_connection

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
        wge_data = query_wge_by_id(event_data['wge_id'])
        grna_class = prepare_guide_rna_class(event_data, wge_data)
        benchling_body = grna_class.as_benchling_req_body(event_data)
        patch_url = benchling_connection.sequence_url + event_data['entity_id'] 
        response = export_to_service(
            benchling_body,
            patch_url,
            benchling_connection.token,
            'patch',
        )

        return response
