from src.rest_calls.send_calls import Caller
from src.utils.exceptions import OligoDirectionInvalid
from src.domain.guideRNA import Oligo
from . import benchling_connection
import json
import sys
sys.path.append("..")


def prepare_oligos_json(oligos, ids):
    return {
        "bases": str(getattr(oligos, 'sequence')),
        # "ids": {
        #     "value": str(getattr(oligos, 'id')),
        # },
        "fields": {

        #     # "bases": {
        #     #     "value": str(getattr(oligos, 'bases')),
        #     # },
        #     # "direction": {
        #         # "value": str(getattr(oligos, 'direction')),
        #     # },
             "Targeton": {
                "value": str(getattr(oligos, 'targeton')),
            },
            "Strand": {
                "value": str(getattr(oligos, 'strand')),
            },
            "Guide RNA": {
                "value": str(getattr(oligos, 'grna'))
            }
        },
        "folderId": str(getattr(oligos, 'folder_id')),
        "name": str(getattr(oligos, 'name')),
        "schemaId": str(getattr(oligos, 'schema_id'))
    }


def export_oligos_to_benchling(oligos):
    benchling_ids = json.load(open('benchling_ids.json'))

    api_caller = Caller(benchling_connection.oligos_url)
    token = benchling_connection.token

    oligos_json = prepare_oligos_json(oligos, benchling_ids)
    print(oligos_json)

    try:
        olgos_id = api_caller.make_request('post', token, oligos_json).json()['id']

    except Exception as err:
        raise Exception(err)

    return olgos_id

def setup_oligo_class(oligo: Oligo, guide_data: dict, benchling_ids: dict, direction: str, name: str = "Guide RNA Oligo", schema_id: str = "ts_wFWXiFSo") -> None:
    oligo.targeton = guide_data["targeton"]
    oligo.folder_id = guide_data["folder_id"]
    oligo.schema_id = schema_id
    oligo.name = name
    if direction == "forward":
        oligo.strand = benchling_ids["forward_strand"]
    elif direction == "reverse":
        oligo.strand = benchling_ids["reverse_strand"]
    else: 
        raise OligoDirectionInvalid(f"Invalid direction given {direction}, expecting \"forward\" or \"reverse\"")
    
    oligo.grna = guide_data["id"]
    return oligo
    
