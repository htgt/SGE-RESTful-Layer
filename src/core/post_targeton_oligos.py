from src.benchling.utils.request_to_benchling import request_to_benchling
import json

def send_targeton_oligo_post_request(body: dict, url: str) -> dict:
    try:
        response = request_to_benchling(
            body,
            url,
            'post'
        )
        return response, 201
    except Exception as err:
        print(json.dumps(err))
        return json.dumps(err), 500
