from src.domain.guideRNA import GuideRNA
from . import benchling_connection, benchling_schema_ids
from src.rest_calls.send_calls import export_to_service_json_response

def patch_guide_rna(guide: GuideRNA, event_data: dict) -> str:
    benchling_body = as_benchling_req_body(guide, event_data)
    patch_url = benchling_connection.sequence_url + '/' + event_data["entity_id"]

    response = export_to_service_json_response(
        benchling_body,
        patch_url,
        benchling_connection.token,
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

    print('Benchling json::', body)
    return body


def get_species_name(ids, species) -> str:
    species_name = "mus_musculus"

    if species == 'Grch37' or species == 'Grch38':
        species_name = "homo_sapience"

    return species_name