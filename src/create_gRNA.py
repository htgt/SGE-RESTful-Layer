from domain.guideRNA import GuideRNA
from rest_calls.send_calls import Caller

folder_id = "lib_BCLzex2p"
grna_schema_id = "ts_61h752hP"
sgrna_schema_id = "ts_3arvRFIl"
positive_strand = "sfso_Mmvf1ki5"
negative_strand = "sfso_zIBEDlMM"

api_caller = Caller()
api_dna_oligos = 'https://tol-sangertest.benchling.com/api/v2/dna-oligos'

data = {
    'id': '1067959134',
    'sequence': 'ACATGGTATTGCAGTAGAC',
    'gene_name': 'A1CF'
}

gRNA = GuideRNA(data)

fwd_sgrna = {
    "bases": str(gRNA.forward_sgRNA()),
    "fields": {
        "Strand": {
            "value": positive_strand,
        }
    },
    "folderId": folder_id,
    "name": f"fwd_{str(getattr(gRNA, 'id'))}",
    "schemaId": sgrna_schema_id
}

fwd_sgrna_json = api_caller.make_post(api_dna_oligos, fwd_sgrna).json()

rev_sgrna = {
    "bases": str(gRNA.reverse_sgRNA()),
    "fields": {
        "Strand": {
            "value": negative_strand,
        }
    },
    "folderId": folder_id,
    "name": f"rev_{str(getattr(gRNA, 'id'))}",
    "schemaId": sgrna_schema_id
}

rev_sgrna_json = api_caller.make_post(api_dna_oligos, rev_sgrna).json()

api_post_data = {
    "bases": str(getattr(gRNA, 'sequence')),
    "fields": {
        "Gene Name": {
            "value": str(getattr(gRNA, 'gene_name')),
        },
        "WGE ID": {
            "value": int(getattr(gRNA, 'id')),
        },
        "Forward sgRNA": {
            "value": str(fwd_sgrna_json['id']),
        },
        "Reverse sgRNA": {
            "value": str(rev_sgrna_json['id']),
        },
    },
    "folderId": folder_id,
    "name": str(getattr(gRNA, 'id')),
    "schemaId": grna_schema_id
}

res = api_caller.make_post(api_dna_oligos, api_post_data)
