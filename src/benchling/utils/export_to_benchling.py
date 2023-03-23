from __future__ import annotations
from typing import TYPE_CHECKING
from src.rest_calls.send_calls import export_to_service, check_response_object
if TYPE_CHECKING:
    from src.benchling import BenchlingConnection

def export_to_benchling(
    json_dict: dict,
    service_url : str,
    connection: BenchlingConnection,
    action : str = 'get',
) -> str:

    response = export_to_service(json_dict, service_url, connection.token, action=action)
    if response.status_code in ["400", "401", "403"] and not response.ok:
        connection.get_store_token()
        response = export_to_service(json_dict, service_url, connection.token, action=action)

    return response

def export_to_benchling_json_response(
    *args,
    **kwargs
) -> str:

    response = export_to_benchling(*args, **kwargs)
    
    json_response = check_response_object(response)

    return json_response