from src.biology.guideRNA import GuideRNA
from . import benchling_connection, benchling_schema_ids
from src.benchling.utils.export_to_benchling import export_to_benchling_json_response

def patch_guide_rna(guide: GuideRNA, event_data: dict) -> str:
    benchling_body = as_benchling_req_body(guide, event_data)
    patch_url = benchling_connection.sequence_url + '/' + event_data["entity_id"]

    response = export_to_benchling_json_response(
        benchling_body,
        patch_url,
        benchling_connection,
        'patch',
    )

    return response['id']

def as_benchling_req_body(guide: GuideRNA, event: dict) -> dict:
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
        'folderId' : event['folder_id'],
        'name': event['name'],
        'schemaId' : event['schema_id'],
    }

    return body
