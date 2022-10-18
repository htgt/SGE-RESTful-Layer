from domain.guideRNA import GuideRNA
import requests
from rest_calls.send_calls import Caller

data = {
    'id': '1067959134',
    'sequence': 'ACATGGTATTGCAGTAGAC',
    'gene_name': 'A1CF'
}

gRNA = GuideRNA(data)

skeleton_post_data = {
    "bases": str(getattr(gRNA, 'sequence')),
    "folderId": "lib_BCLzex2p",
    "name": str(getattr(gRNA, 'id'))
}

api_post_data = {
    "bases": str(getattr(gRNA, 'sequence')),
    "fields": {
        #"gene_name": {
        #    "value": str(getattr(gRNA, 'gene_name')),
        #},
        "wge_id": {
            "value": int(getattr(gRNA, 'id')),
        },
        #"forward_sgrna": {
        #    "value": str((gRNA.forward_sgRNA())),
        #},
        #"reverse_sgrna": {
        #    "value": str(gRNA.reverse_sgRNA()),
        #},
    },
    "folderId": "lib_BCLzex2p",
    "name": str(getattr(gRNA, 'id')),
    "schemaId": "ts_61h752hP"
}

api_caller = Caller()

res = api_caller.make_post('https://tol-sangertest.benchling.com/api/v2/dna-oligos', api_post_data)

if res.ok:
    print(res.status_code)
else:
    print(res.status_code)
    print(res.text)
