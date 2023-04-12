from src.rest_calls.send_calls import Caller
import json
import sys
sys.path.append("..")


def prepare_sgrna_json(gRNA, strand, ids):
    if strand == '+':
        bases = str(gRNA.forward_sgRNA())
        strand_id = ids['positive_strand']
        name = f"fwd_{str(getattr(gRNA, 'id'))}"
    elif strand == '-':
        bases = str(gRNA.reverse_sgRNA())
        strand_id = ids['negative_strand']
        name = f"rev_{str(getattr(gRNA, 'id'))}"
    else:
        return None

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


def export_grna_to_benchling(gRNA, benchling_connection):
    benchling_ids = json.load(open('benchling_ids.json'))

    api_caller = Caller(benchling_connection.oligos_url)
    token = benchling_connection.token

    fwd_sgrna = prepare_sgrna_json(gRNA, '+', benchling_ids)
    rev_sgrna = prepare_sgrna_json(gRNA, '-', benchling_ids)

    try:
        fwd_sgrna_id = api_caller.make_request('post', token, fwd_sgrna).json()['id']
        rev_sgrna_id = api_caller.make_request('post', token, rev_sgrna).json()['id']

        api_post_data = prepare_grna_json(gRNA, fwd_sgrna_id, rev_sgrna_id, benchling_ids)

        grna_id = api_caller.make_request('post', token, api_post_data).json()['id']

    except Exception as err:
        raise Exception(err)

    return [fwd_sgrna_id, rev_sgrna_id]
