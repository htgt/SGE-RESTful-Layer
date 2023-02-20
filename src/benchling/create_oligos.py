from src.rest_calls.send_calls import Caller
from . import benchling_connection
import json
import sys
sys.path.append("..")


def prepare_oligos_json(gRNA, fwd_sgrna_id, rev_sgrna_id, ids):
    return {
        "bases": str(getattr(gRNA, 'sequence')),
        "fields": {
            "Gene Name": {
                "value": str(getattr(gRNA, 'gene_name')),
            },
            "WGE ID": {
                "value": int(getattr(gRNA, 'id')),
            },
            "Forward sgRNA": {
                "value": str(fwd_sgrna_id),
            },
            "Reverse sgRNA": {
                "value": str(rev_sgrna_id),
            },
        },
        "folderId": ids['folder_id'],
        "name": str(getattr(gRNA, 'id')),
        "schemaId": ids['grna_schema_id']
    }


def export_oligos_to_benchling(oligos):
    benchling_ids = json.load(open('benchling_ids.json'))

    api_caller = Caller(benchling_connection.oligos_url)
    token = benchling_connection.token

    oligos_json = prepare_oligos_json(oligos, '+', benchling_ids)

    try:
        fwd_sgrna_id = api_caller.make_request('post', token, oligos_json).json()['id']

    except Exception as err:
        raise Exception(err)

    return fwd_sgrna_id
