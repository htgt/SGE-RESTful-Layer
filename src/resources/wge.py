from flask import request
from flask_restful import Resource

#from src.wge.wge import transform_wge_event, patch_wge_data_to_service

import requests
import json

class WGEEndpoint(Resource):
    def get(self):
       # wge_id = request.args.get('id')
       # wge_json = query_wge_by_id(wge_id)

       # return wge_json
        return 'WGE', 200
    
    def post(self):
        data = request.json
    #    event_data = transform_wge_event(data)
    #    response = patch_wge_data_to_service(event_data)

        return data, 200
