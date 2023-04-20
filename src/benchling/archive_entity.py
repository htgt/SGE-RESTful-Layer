from src.benchling.utils.request_to_benchling import request_to_benchling


def archive_oligo(id, url):
    return archive_entity(id, 'dnaOligo', url)

def archive_entity(entity_id, entity_type, url):
    json = prepare_archive_json(entity_type, entity_id)

    result = send_archive_request(url, json)

    return result

def send_archive_request(url, json):
    response = request_to_benchling(url, 'post', json=json)

    return response

def prepare_archive_json(entity_type, entity_id):
    return {
        generate_entity_ids_str(entity_type): [entity_id],
        "reason": "Made in error"
    }

def generate_entity_ids_str(entity_type):
    return entity_type + "Ids"