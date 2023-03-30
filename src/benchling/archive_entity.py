from src.rest_calls.send_calls import Caller
from . import benchling_connection


def archive_oligo(id):
    return archive_entity(id, 'dnaOligo', benchling_connection.oligos_url)

def archive_entity(entity_id, entity_type, url):
    json = prepare_archive_json(entity_type, entity_id)

    result = send_archive_request(url, json)

    return result

def send_archive_request(url, json):
    api_caller = Caller(url)
    token = benchling_connection.token

    response = api_caller.make_request('post', token, json)

    return response

def prepare_archive_json(entity_type, entity_id):
    return {
        generate_entity_ids_str(entity_type): [entity_id],
        "reason": "Made in error"
    }

def generate_entity_ids_str(entity_type):
    return entity_type + "Ids"