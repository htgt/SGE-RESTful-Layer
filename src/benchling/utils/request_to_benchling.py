from src.rest_calls.send_calls import request_to_service, check_response_object
import requests

def request_to_benchling(
    service_url : str,
    action : str ,
    json = {}
) -> requests.Response:
    from src.benchling.connection.benchling_connection import benchling_connection
    response = request_to_service(json, service_url, benchling_connection.token, action=action)
    if response.status_code in [400, 401, 403] and not response.ok:
        print("Regenerating token...")
        benchling_connection.get_store_token()
        response = request_to_service(json, service_url, benchling_connection.token, action=action)

    return response


def request_to_benchling_json_response(
    *args,
    **kwargs
) -> dict:

    response = request_to_benchling(*args, **kwargs)

    json_response = check_response_object(response)

    return json_response
