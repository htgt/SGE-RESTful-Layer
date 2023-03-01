from flask import request
from flask_restful import Resource
from src.wge.wge import query_wge_by_id
import requests


class WGEEndpoint(Resource):
    def __transform_event_input_data(data):
        data_entity = data['detail']['entity']
        
        wge_grna_data = {}
        wge_grna_data['entity_id'] = data_entity['id']
        wge_grna_data['wge_id'] = data_entity['fields']['WGE_ID']['value']
        wge_grna_data['targeton_id'] = data_entity['fields']['Targeton']['value']

        return wge_grna_data

    def get(self):
        wge_id = request.args.get('id')
        wge_json = query_wge_by_id(wge_id)
        return wge_json
    
    def post(self, id):
        data = request.json
        print(data)
        return data
        #return (export_return_forward, export_return_reverse), 200

            #except Exception as err:
            #    return json.dumps(err), 500
       # else:
       #     return "Incorrect input data", 404
