from __future__ import annotations
from dataclasses import dataclass
from typing import List, Tuple
from src.utils.base_classes import BaseClass

from src.utils.exceptions import OligoDirectionInvalid
from src.biology.guideRNA import Oligo, GuideRNAOligos
from src.benchling import benchling_schema_ids
import sys
from src.utils.base_classes import BaseClass
from src.benchling.utils.request_to_benchling import request_to_benchling_json_response
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
class BenchlingOligosPair(BaseClass):
    forward: BenchlingOligo
    reverse: BenchlingOligo

    def to_list_dicts(self) -> List[dict]:
        list_of_dicts = []
        for field in self.get_fields():
            return_dict = getattr(self, field)._asdict()
            list_of_dicts.append(return_dict)
        return list_of_dicts


def prepare_oligo_json(oligo: BenchlingOligo) -> dict:
    return {
        "bases": str(getattr(oligo, 'sequence')),
        "fields": {
            "Targeton": {
                "value": str(getattr(oligo, 'targeton')),
            },
            "Strand": {
                "value": str(getattr(oligo, 'strand')),
            },
            "Guide RNA": {
                "value": str(getattr(oligo, 'grna'))
            }
        },
        "isCircular": False,
        "folderId": str(getattr(oligo, 'folder_id')),
        "name": str(getattr(oligo, 'name')),
        "schemaId": str(getattr(oligo, 'schema_id'))
    }


def export_oligos_to_benchling(oligos: BenchlingOligosPair, url: str) -> Tuple[str, str]:
    oligo_forward_json = prepare_oligo_json(oligos.forward)
    oligo_reverse_json = prepare_oligo_json(oligos.reverse)

    oligo_forward = request_to_benchling_json_response(
        oligo_forward_json,
        url,
        'post',
    )
    oligo_reverse = request_to_benchling_json_response(
        oligo_reverse_json,
        url,
        'post',
    )


    return (oligo_forward, oligo_reverse)


def setup_oligo_pair_class(oligos: GuideRNAOligos, guide_data: dict) -> BenchlingOligosPair:
    benchling_ids = benchling_schema_ids.ids
    forward = setup_oligo_class(
        oligos.forward,
        guide_data,
        benchling_ids,
        'forward',
    )
    reverse = setup_oligo_class(
        oligos.reverse,
        guide_data,
        benchling_ids,
        'reverse',
    )
    benchling_oligo_pair = BenchlingOligosPair(forward, reverse)
    return benchling_oligo_pair

def setup_oligo_class(
        oligo: Oligo,
        guide_data: dict,
        benchling_ids: dict,
        direction: str,
        name: str = "Guide RNA Oligo"
) -> BenchlingOligo:
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
