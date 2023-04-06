from src.benchling.utils.export_to_benchling import export_to_benchling
from src.benchling.connection.benchling_connection import benchling_connection
import json

def send_targeton_oligo_post_request(body: dict) -> dict:
    try:
        response = export_to_benchling(
            body,
            benchling_connection.custom_entity_url,
            benchling_connection,
            'post'
        )
        return response, 201
    except Exception as err:
        print(json.dumps(err))
        return json.dumps(err), 500
