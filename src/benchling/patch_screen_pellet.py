from src.benchling import benchling_schema_ids, benchling_urls
from src.benchling.utils.request_to_benchling import request_to_benchling_json_response


def patch_screen_pellet(pellet_data: dict, pellet_id: str) -> dict:
    url = benchling_urls.custom_entity_url
    benchling_body = as_benchling_req_body(pellet_data)
    patch_url = url + '/' + pellet_id

    json_response = request_to_benchling_json_response(
        patch_url,
        'patch',
        benchling_body
    )

    return json_response


def as_benchling_req_body(pellet_data: dict) -> dict:
    schema_id = benchling_schema_ids.ids['schemas']['screen_pellet_schema_id']
    status = convert_status_id_to_status(pellet_data['run_status'])
    status_dropdown_id = benchling_schema_ids.ids['dropdowns']['sequencing_status'][status]
    run_id = pellet_data['id_run']

    body = {
        'fields': {
            'iRODs Path': {
                'value': pellet_data['irods_data_relative_path'],
            },
            'Sequencing Status': {
                'value': status_dropdown_id,
            },
            'Run ID': {
                'value': run_id,
            },
        },
        'schemaId': schema_id,
    }

    return body


def convert_status_id_to_status(status_id: int) -> str:
    if status_id == 20:
        status = 'finished'
    else:
        status = 'submitted'
    return status
