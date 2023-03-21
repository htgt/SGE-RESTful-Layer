from dataclasses import dataclass
from typing import List
from src.utils.base_classes import BaseClass

from src.rest_calls.send_calls import Caller, export_to_service_json_response
from src.utils.exceptions import OligoDirectionInvalid
from src.biology.guideRNA import Oligo
from . import BenchlingConnection, benchling_schema_ids
import json
import sys

sys.path.append("..")


@dataclass
class BenchlingOligo(Oligo, BaseClass):
    targeton: str
    folder_id: str
    schema_id: str
    name: str
    strand: str
    grna: str


@dataclass
class OligosPair(BaseClass):
    forward: BenchlingOligo
    reverse: BenchlingOligo

    def to_list_dicts(self) -> List[dict]:
        list_of_dicts = []
        for field in self.get_fields():
            return_dict = getattr(self, field)._asdict()
            list_of_dicts.append(return_dict)
        return list_of_dicts


def prepare_oligo_json(oligos: OligosPair) -> dict:
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
        "isCircular": False,
        "folderId": str(getattr(oligos, 'folder_id')),
        "name": str(getattr(oligos, 'name')),
        "schemaId": str(getattr(oligos, 'schema_id'))
    }


def export_oligos_to_benchling(oligos: OligosPair, benchling_connection: BenchlingConnection):
    oligo_forward_json = prepare_oligo_json(oligos.forward)
    oligo_reverse_json = prepare_oligo_json(oligos.reverse)

    oligo_forward = export_to_service_json_response(
        oligo_forward_json,
        benchling_connection.sequence_url,
        benchling_connection.token,
        'post',
    )
    oligo_reverse = export_to_service_json_response(
        oligo_reverse_json,
        benchling_connection.sequence_url,
        benchling_connection.token,
        'post',
    )

    return (oligo_forward['id'], oligo_reverse['id'])


def setup_oligo_pair_class(oligos: OligosPair, guide_data: dict) -> OligosPair:
    benchling_ids = benchling_schema_ids.ids

    oligos.forward = setup_oligo_class(
        oligos.forward,
        guide_data,
        benchling_ids,
        'forward',
    )
    oligos.reverse = setup_oligo_class(
        oligos.reverse,
        guide_data,
        benchling_ids,
        'reverse',
    )
    return oligos

def setup_oligo_class(
        oligo: Oligo,
        guide_data: dict,
        benchling_ids: dict,
        direction: str,
        name: str = "Guide RNA Oligo"
) -> None:
    if direction == "forward":
        strand = benchling_ids["dropdowns"]["sense"]
    elif direction == "reverse":
        strand = benchling_ids["dropdowns"]["antisense"]
    else:
        raise OligoDirectionInvalid(
            f"Invalid direction given {direction}, expecting \"forward\" or \"reverse\"")

    benchling_oligo = BenchlingOligo(
        sequence=oligo.sequence,
        targeton=guide_data["targeton"],
        folder_id=guide_data["folder_id"],
        schema_id=benchling_ids["schemas"]["grna_oligo_schema_id"],
        name=name,
        strand=strand,
        grna=guide_data["id"]
    )

    return benchling_oligo
