import json
from src.domain.guideRNA import GuideRNAOligo
from src.benchling import benchling_connection
from src.benchling.create_oligos import export_oligos_to_benchling, setup_oligo_pair_class
from src.benchling.get_sequence import get_sequence
from src.wge.wge import patch_wge_data_to_service, transform_wge_event

BENCHLING_GUIDE_RNA_SCHEMA_ID = "ts_vGZYroiQ"
BENCHLING_ENTITY_REGISTERED_EVENT = "v2.entity.registered"

def handle_guide_event(data : dict) -> dict:
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
    guide_data = transform_event_input_data(data)
    guide_data["seq"] = get_sequence(guide_data["id"])
    oligos = GuideRNAOligo(guide_data["seq"]).create_oligos()
    benchling_ids = json.load(open('benchling_ids.json'))
    oligos = setup_oligo_pair_class(oligos, guide_data, benchling_ids)
    
    return oligos

def check_wge_id(data : dict) -> bool:
    check = False
    if data['detail']['entity']['fields']['WGE ID']:
        check = True
    return check

def check_event_is_guide_rna(data: dict) -> bool:
    bool_check = True
    if not data["detail-type"] == BENCHLING_ENTITY_REGISTERED_EVENT:
        bool_check = False
    if not data["detail"]["entity"]["schema"]["id"] == BENCHLING_GUIDE_RNA_SCHEMA_ID:
        bool_check = False

def transform_event_input_data(data):
    guide_data = {}

    guide_data["id"] = data["detail"]["entity"]["id"]
    guide_data["targeton"] = data["detail"]["entity"]["fields"]["Targeton"]["value"]
    guide_data["folder_id"] = data["detail"]["entity"]["folderId"]
    guide_data["schemaid"] = "ts_wFWXiFSo"
    guide_data["name"] = "Guide RNA Oligo"

    return guide_data
