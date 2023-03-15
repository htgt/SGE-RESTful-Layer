import json
from src.domain.guideRNA import GuideRNAOligo
from src.benchling import benchling_connection, benchling_schema_ids

from src.benchling.create_oligos import export_oligos_to_benchling, setup_oligo_pair_class
from src.benchling.get_sequence import get_sequence
from src.wge.wge import patch_wge_data_to_service, transform_wge_event


def handle_guide_event(data : dict) -> dict:
    response = {}
    try:
        response['grna'] = patch_grna_event(data)
    except Exception as err:
        return json.dumps(err), 500
    try:
        response['oligos'] = post_grna_oligos_event(data)
        return response, 201
    except Exception as err: 
        return json.dumps(err), 500

def patch_grna_event(data : dict) -> dict:
    if check_wge_id(data):
        wge_event = transform_wge_event(data)
        response = patch_wge_data_to_service(wge_event)

        return response

def post_grna_oligos_event(data : dict) -> dict:
    if check_event_is_guide_rna(data):
        oligos = transform_grna_oligos(data)
        export_response = export_oligos_to_benchling(
            oligos,
            benchling_connection
        )
        return export_response
    else:
        return "Incorrect input data", 404

def transform_grna_oligos(data : dict) -> dict:
    benchling_ids = benchling_schema_ids.ids

    guide_data = transform_event_input_data(data, benchling_ids)
    guide_data["seq"] = get_sequence(guide_data["id"])
    oligos = GuideRNAOligo(guide_data["seq"]).create_oligos()
    #benchling_ids = json.load(open('benchling_ids.json'))

    oligos = setup_oligo_pair_class(oligos, guide_data, benchling_ids)
    
    return oligos

def check_wge_id(data : dict) -> bool:
    check = False
    if data['detail']['entity']['fields']['WGE ID']:
        check = True
    return check

def check_event_is_guide_rna(data: dict) -> bool:
    bool_check = True
    if not data["detail-type"] == benchling_schema_ids.ids["events"]["entity_registered"]:
        bool_check = False
    if not data["detail"]["entity"]["schema"]["id"] == benchling_schema_ids.ids["schemas"]["grna_schema_id"]:
        bool_check = False
    return bool_check


def transform_event_input_data(data, ids):
    guide_data = {}

    guide_data["id"] = data["detail"]["entity"]["id"]
    guide_data["targeton"] = data["detail"]["entity"]["fields"]["Targeton"]["value"]
    guide_data["folder_id"] = data["detail"]["entity"]["folderId"]
    guide_data["schemaid"] = ids["schemas"]["grna_oligo_schema_id"]
    guide_data["name"] = "Guide RNA Oligo"

    return guide_data
