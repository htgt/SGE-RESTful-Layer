from __future__ import annotations
from typing import TYPE_CHECKING
from src.biology.guideRNA import GuideRNA
from src.benchling.utils.export_to_benchling import export_to_benchling_json_response

if TYPE_CHECKING:
    from src.benchling.connection.connection_class import BenchlingConnection
    from src.benchling import BenchlingSchemaIds

def patch_guide_rna(guide: GuideRNA, event_data: dict, benchling_connection: BenchlingConnection, benchling_schema_ids: BenchlingSchemaIds) -> str:
    benchling_body = as_benchling_req_body(guide, event_data, benchling_schema_ids)
    patch_url = benchling_connection.sequence_url + '/' + event_data["entity_id"]

    response = export_to_benchling_json_response(
        benchling_body,
        patch_url,
        benchling_connection,
        'patch',
    )

    return response

def as_benchling_req_body(guide: GuideRNA, event: dict, benchling_schema_ids: BenchlingSchemaIds) -> dict:
    species_benchling_id = benchling_schema_ids.ids["dropdowns"]["species"][guide.species]

    body = {
        'bases': guide.sequence,
        'fields'   : {
            'WGE ID'       : {
                'value': guide.wge_id,
            },
            'Targeton'               : {
                'value': event['targeton_id'],
            },
            #    'Strand' : {
            #        'value' : self.strand,
            #    },
            'WGE Hyperlink'          : {'value': guide.wge_link, },
            'Off Target Summary Data': {'value': guide.off_targets, },
            'Species'                : {'value': species_benchling_id, },
        },
        'name': event['name'],
        'schemaId' : event['schema_id'],
    }

    return body
