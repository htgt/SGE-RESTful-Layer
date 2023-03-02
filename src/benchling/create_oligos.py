from src.rest_calls.send_calls import Caller
from src.utils.exceptions import OligoDirectionInvalid
from src.domain.guideRNA import Oligo
from dataclasses import dataclass
from . import BenchlingConnection
import json
import sys
from src.utils.classes import BaseClass
sys.path.append("..")


@dataclass
class BenchlingOligo(Oligo, BaseClass):
    targeton: str
    folder_id: str
    schema_id: str
    name: str
    strand: str
    grna: str


def prepare_oligos_json(oligos, ids):
    return {
        "bases": str(getattr(oligos, 'sequence')),
        "fields": {
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


def export_oligos_to_benchling(oligos: BenchlingOligo, benchling_connection: BenchlingConnection, benchling_ids_path='benchling_ids.json'):
    benchling_ids = json.load(open(benchling_ids_path))

    api_caller = Caller(benchling_connection.oligos_url)
    token = benchling_connection.token

    oligos_json = prepare_oligos_json(oligos, benchling_ids)

    try:
        olgos_id = api_caller.make_request('post', token, oligos_json).json()['id']

    except Exception as err:
        raise Exception(err)

    return olgos_id


def setup_oligo_class(oligo: Oligo, guide_data: dict, benchling_ids: dict, direction: str, name: str = "Guide RNA Oligo", schema_id: str = "ts_wFWXiFSo") -> None:
    if direction == "forward":
        strand = benchling_ids["forward_strand"]
    elif direction == "reverse":
        strand = benchling_ids["reverse_strand"]
    else:
        raise OligoDirectionInvalid(
            f"Invalid direction given {direction}, expecting \"forward\" or \"reverse\"")

    benchling_oligo = BenchlingOligo(
        sequence=oligo.sequence,
        targeton=guide_data["targeton"],
        folder_id=guide_data["folder_id"],
        schema_id=schema_id,
        name=name,
        strand=strand,
        grna=guide_data["id"]
    )

    return benchling_oligo
