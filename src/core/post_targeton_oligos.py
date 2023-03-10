from src.rest_calls.send_calls import export_to_service
from src.benchling import benchling_connection

def send_targeton_oligo_post_request(body: dict) -> dict:
    try:
        response = export_to_service(
            body,
            benchling_connection.custom_entity_url,
            benchling_connection,
            'post'
        )
        return response, 201
    except Exception as err:
        print(json.dumps(err))
        return json.dumps(err), 500
