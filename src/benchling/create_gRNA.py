import json

from benchling.auth_uitls import AuthUtils
from rest_calls.send_calls import Caller
from domain.guideRNA import GuideRNA

def prepare_sgrna_json(gRNA, strand, ids):
    if strand == '+':
        bases = str(gRNA.forward_sgRNA())
        strand_id = ids['positive_strand']
        name = f"fwd_{str(getattr(gRNA, 'id'))}"
    else:
        bases = str(gRNA.reverse_sgRNA())
        strand_id = ids['negative_strand']
        name = f"rev_{str(getattr(gRNA, 'id'))}"

    return {
        "bases": bases,
        "fields": {
            "Strand": {
                "value": strand_id,
            }
        },
        "folderId": ids['folder_id'],
        "name": name,
        "schemaId": ids['sgrna_schema_id']
    }

def prepare_grna_json(gRNA, fwd_sgrna_id, rev_sgrna_id, ids):
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

def export_grna_to_benchling(data):
    benchling_ids = json.load(open('benchling_ids.json'))

    gRNA = GuideRNA(data)
    api_caller = Caller('https://tol-sangertest.benchling.com/api/v2/dna-oligos')

    fwd_sgrna = prepare_sgrna_json(gRNA, '+', benchling_ids)
    rev_sgrna = prepare_sgrna_json(gRNA, '-', benchling_ids)

    fwd_sgrna_id = api_caller.make_post(AuthUtils.get_access_token, fwd_sgrna).json()['id']
    rev_sgrna_id = api_caller.make_post(AuthUtils.get_access_token, rev_sgrna).json()['id']

    api_post_data = prepare_grna_json(gRNA, fwd_sgrna_id, rev_sgrna_id, benchling_ids)

    api_caller.make_post(AuthUtils.get_access_token, api_post_data)
