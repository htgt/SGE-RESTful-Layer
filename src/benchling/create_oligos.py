from src.rest_calls.send_calls import Caller
from src.utils.exceptions import OligoDirectionInvalid
from src.domain.guideRNA import Oligo, OligosPair
from dataclasses import dataclass
from . import BenchlingConnection
import json
import sys
from src.utils.base_classes import BaseClass
from src.benchling.auth_utils import export_to_benchling
sys.path.append("..")


@dataclass
class BenchlingOligo(Oligo, BaseClass):
    targeton: str
    folder_id: str
    schema_id: str
    name: str
    strand: str
    grna: str


def prepare_oligo_json(oligos: BenchlingOligo) -> dict:
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


def export_oligos_to_benchling(oligos: BenchlingOligo, benchling_connection: BenchlingConnection):
    oligo_forward_json = prepare_oligo_json(oligos.forward)
    oligo_reverse_json = prepare_oligo_json(oligos.reverse)
    oligo_forward = export_to_benchling(
        oligo_forward_json,
        benchling_connection.oligos_url,
        benchling_connection,
        'post',
    )
    oligo_reverse = export_to_benchling(
        oligo_reverse_json,
        benchling_connection.oligos_url,
        benchling_connection,
        'post',
    )

    return (oligo_forward['id'], oligo_reverse['id'])


def setup_oligo_pair_class(oligos: OligosPair, guide_data: dict, benchling_ids: dict) -> OligosPair:
    # Foward
    oligos.forward = setup_oligo_class(
        oligos.forward,
        guide_data,
        benchling_ids,
        'forward',
    )
    # Reverse
    oligos.reverse = setup_oligo_class(
        oligos.reverse,
        guide_data,
        benchling_ids,
        'reverse',
    )
    return oligos


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
