import json

from . import benchling_connection, benchling_schema_ids
from src.benchling.utils.export_to_benchling import export_to_benchling_json_response


def patch_screen_pellet(pellet_data: dict, pellet_id: str) -> dict:
    benchling_body = as_benchling_req_body(pellet_data)
    patch_url = benchling_connection.custom_entity_url + '/' + pellet_id

    json_response = export_to_benchling_json_response(
        benchling_body,
        patch_url,
        benchling_connection,
        'patch',
    )

    return json_response

def as_benchling_req_body(pellet_data: dict) -> dict:
    schema_id = benchling_schema_ids.ids['schemas']['screen_pellet_schema_id']
    status = benchling_schema_ids.ids['dropdowns']['sequencing_status'][pellet_data['run_status']]

    body = {
        'fields': {
            'iRODs Path': {
                'value': pellet_data['irods_data_relative_path'],
            },
            'Sequencing Status': {
                'value': status,
            },
        },
        'schemaId' : schema_id,
    }

    return body