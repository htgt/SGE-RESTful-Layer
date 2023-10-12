import json
from src.biology.guideRNA import GuideRNAOligos
from src.benchling.create_oligos import export_oligos_to_benchling, setup_oligo_pair_class
from src.benchling.get_sequence import get_sequence
from src.benchling.patch_guide_rna import patch_guide_rna
from src.wge.wge import query_wge_by_id, transform_wge_event, prepare_guide_rna_class
from src.benchling import benchling_schema_ids



def handle_guide_event(data : dict) -> dict:
    response = {}
    try:
        print(data['detail']['entity']['fields']['WGE ID']['value'])
        wge_response = query_wge_by_id(data['detail']['entity']['fields']['WGE ID']['value'])
        print('WGE response', wge_response)

        response['grna'] = patch_grna_event(data, wge_response)
        print('GRNA patched')
    except Exception as err:
        return str(err), 500
    try:
        print('Try POST grna oligos')

        response['oligos'] = post_grna_oligos_event(data)
        return response, 201
    except Exception as err: 
        return str(err), 500


def patch_grna_event(data : dict) -> dict:
    wge_event = transform_wge_event(data)
    response = patch_wge_data_to_service(wge_event)

    return response

def patch_wge_data_to_service(event_data : dict, wge_data: dict) -> dict:
    grna_object = prepare_guide_rna_class(event_data, wge_data)

    response = patch_guide_rna(grna_object, event_data)

    return response


def post_grna_oligos_event(data : dict) -> dict:
    print('Start post_grna_oligos_event', data)
    oligos = transform_grna_oligos(data)

    print('Calculated Oligos:', oligos)

    export_response = export_oligos_to_benchling(oligos)

    return export_response



def transform_grna_oligos(data : dict) -> dict:
    print('Start transform_grna_oligos')

    benchling_ids = benchling_schema_ids.ids

    guide_data = transform_event_input_data(data, benchling_ids)

    print('Get Sequence', guide_data)

    guide_data["seq"] = get_sequence(guide_data["id"])
    oligos = GuideRNAOligos(guide_data["seq"])

    oligos = setup_oligo_pair_class(oligos, guide_data)
    
    return oligos


def check_wge_id(data : dict) -> bool:
    check = False
    if data['detail']['entity']['fields']['WGE ID']:
        check = True
    return check


def transform_event_input_data(data: dict, ids: dict) -> dict:
    guide_data = {}

    guide_data["id"] = data["detail"]["entity"]["id"]
    guide_data["targeton"] = data["detail"]["entity"]["fields"]["Targeton"]["value"]
    guide_data["folder_id"] = data["detail"]["entity"]["folderId"]
    guide_data["schemaid"] = ids["schemas"]["grna_oligo_schema_id"]
    guide_data["name"] = "Guide RNA Oligo"

    return guide_data
