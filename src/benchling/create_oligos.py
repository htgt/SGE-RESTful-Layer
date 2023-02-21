from src.rest_calls.send_calls import Caller
from . import benchling_connection
import json
import sys
sys.path.append("..")

def prepare_oligos_json(oligos, ids):
    return {
        "bases": str(getattr(oligos, 'sequence')),
        "fields": {
            "Oligos": {
                "value": str(getattr(oligos, 'gene_name')),
            },
            "ID": {
                "value": int(getattr(oligos, 'id')),
            },
            "Targeton": {
                "value": str(getattr(oligos, 'targeton')),
            },
            "Sequence": {
                "value": str(getattr(oligos, 'sequence')),
            },
        },
        "folderId": ids['folder_id'],
        "name": str(getattr(oligos, 'id')),
        "schemaId": ids['grna_schema_id']
    }


def export_oligos_to_benchling(oligos):
    benchling_ids = json.load(open('benchling_ids.json'))

    api_caller = Caller(benchling_connection.oligos_url)
    token = benchling_connection.token

    oligos_json = prepare_oligos_json(oligos, '+', benchling_ids)

    try:
        olgos_id = api_caller.make_request('post', token, oligos_json).json()['id']

    except Exception as err:
        raise Exception(err)

    return olgos_id
