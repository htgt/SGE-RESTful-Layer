import json

from . import benchling_connection, benchling_schema_ids
from src.benchling.utils.export_to_benchling import export_to_benchling_json_response


def patch_screen_pellet(pellet_data: json, pellet_id: str) -> str:
    benchling_body = as_benchling_req_body(pellet_data)
    patch_url = benchling_connection.custom_entity_url + '/' + pellet_id

    response = export_to_benchling_json_response(
        benchling_body,
        patch_url,
        benchling_connection,
        'patch',
    )

    return response

def as_benchling_req_body(pellet_data: json) -> dict:
    SCHEMA_ID = 'ts_agkeO1Jd'

    body = {
        'fields': {
            'iRODs Path': {
                'value': pellet_data['irods_data_relative_path'],
            },
            'Sequencing Status': {
                'value': pellet_data['run_status'],
            },
        },
        'schemaId' : SCHEMA_ID,
    }

    return body