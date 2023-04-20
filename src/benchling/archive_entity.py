from src.rest_calls.send_calls import Caller


def archive_oligo(id, benchling_connection):
    return archive_entity(id, 'dnaOligo', benchling_connection.oligos_url, benchling_connection.token)

def archive_entity(entity_id, entity_type, url, token):
    json = prepare_archive_json(entity_type, entity_id)

    result = send_archive_request(url, json, token)

    return result

def send_archive_request(url, json, token):
    api_caller = Caller(url)

    response = api_caller.make_request('post', token, json)

    return response

def prepare_archive_json(entity_type, entity_id):
    return {
        generate_entity_ids_str(entity_type): [entity_id],
        "reason": "Made in error"
    }

def generate_entity_ids_str(entity_type):
    return entity_type + "Ids"