from flask import request
from flask_restful import Resource
from src.domain.guideRNA import GuideRNAOligo
from src.benchling.get_sequence import get_sequence
from src.benchling.create_oligos import export_oligos_to_benchling, setup_oligo_pair_class
from src.benchling import benchling_connection
import json

BENCHLING_GUIDE_RNA_SCHEMA_ID = "ts_vGZYroiQ"
BENCHLING_ENTITY_REGISTERED_EVENT = "v2.entity.registered"


class GuideEndpoint(Resource):
    def __transform_event_input_data(self, data):
        guide_data = {}

        guide_data["id"] = data["detail"]["entity"]["id"]
        guide_data["targeton"] = data["detail"]["entity"]["fields"]["Targeton"]["value"]
        guide_data["folder_id"] = data["detail"]["entity"]["folderId"]
        guide_data["schemaid"] = "ts_wFWXiFSo"
        guide_data["name"] = "Guide RNA Oligo"

        return guide_data

    def get(self, id):

        return id, 201

    def post(self):
        data = request.json

        if check_event_is_guide_rna(data):
            try:
                guide_data = self.__transform_event_input_data(data)
                guide_data["seq"] = get_sequence(guide_data["id"])
                oligos = GuideRNAOligo(guide_data["seq"]).create_oligos()
                benchling_ids = json.load(open('benchling_ids.json'))
                oligos = setup_oligo_pair_class(oligos, guide_data, benchling_ids)

                export_return = export_oligos_to_benchling(
                    oligos, 
                    benchling_connection
                )

                return export_return, 200

            except Exception as err:
                return json.dumps(err), 500
        else:
            return "Incorrect input data", 404


def check_event_is_guide_rna(data: dict) -> bool:
    bool_check = True
    if not data["detail-type"] == BENCHLING_ENTITY_REGISTERED_EVENT:
        bool_check = False
    if not data["detail"]["entity"]["schema"]["id"] == BENCHLING_GUIDE_RNA_SCHEMA_ID:
        bool_check = False

    return bool_check
