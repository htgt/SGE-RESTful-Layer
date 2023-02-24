from src.rest_calls.send_calls import Caller
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
