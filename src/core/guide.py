from src.resources.wge import transform_wge_event
from src.wge.wge import patch_wge_data_to_service

def patch_grna_event(data : dict) -> dict:
    if check_wge_id(data):
        wge_event = transform_wge_event(data)
        response = patch_wge_data_to_service(wge_event)

        return response

def check_wge_id(data : dict) -> bool:
    check = False
    if data['detail']['entity']['fields']['WGE ID']:
        check = True
    return check
