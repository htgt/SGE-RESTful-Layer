from flask import request
from flask_restful import Resource
from src.domain.guideRNA import GuideRNAOligo
from src.benchling.get_sequence import get_sequence
from src.benchling.create_oligos import prepare_oligos_json, export_oligos_to_benchling
import json

BENCHLING_GUIDE_RNA_SCHEMA_ID = "ts_vGZYroiQ"
BENCHLING_ENTITY_REGISTERED_EVENT = "v2.entity.registered"

class GuideEndpoint(Resource):
    def __transform_event_input_data(self,data):
        guide_data = {}

        guide_data["id"] = data["detail"]["entity"]["id"]
        # guide_data["seq"] = data["detail"]["entity"]["fields"]["Guide Sequence"]["value"]
        guide_data["targeton"] = data["detail"]["entity"]["fields"]["Targeton"]["value"]
        guide_data["folder_id"] = data["detail"]["entity"]["folderId"]

        return guide_data
    
    def get(self, id):
        
        return id, 201
    
    def post(self):
        data = request.json

        if check_event_is_guide_rna(data):
            try:
                guide_data = self.__transform_event_input_data(data)
                guide_data["seq"] = get_sequence(guide_data["id"])
                # print(guide_data)
                oligos = GuideRNAOligo(guide_data["seq"]).create_oligos()
                print(oligos)
                
                oligos.forward.targeton = guide_data["targeton"]
                oligos.forward.folder_id = guide_data["folder_id"]
                # oligos.forward.id = guide_data["id"]
                oligos.forward.schema_id = "ts_wFWXiFSo"
                oligos.forward.name = "Guide RNA Oligo"
                oligos.reverse.targeton = guide_data["targeton"]
                oligos.reverse.folder_id = guide_data["folder_id"]
                # oligos.reverse.id = guide_data["id"]
                oligos.reverse.schema_id = "ts_wFWXiFSo"
                oligos.reverse.name = "Guide RNA Oligo"
                export_return = export_oligos_to_benchling(oligos.forward)
                export_return = export_oligos_to_benchling(oligos.reverse)
                
                return export_return, 200

            except Exception as err:
                return json.dumps(err), 500
        else:
            return "Incorrect input data", 404
    
def check_event_is_guide_rna(data:dict) -> bool:
    bool_check = True
    if not data["detail-type"] == BENCHLING_ENTITY_REGISTERED_EVENT:
        bool_check = False
    if not data["detail"]["entity"]["schema"]["id"] == BENCHLING_GUIDE_RNA_SCHEMA_ID:
        bool_check = False
        
    return bool_check
