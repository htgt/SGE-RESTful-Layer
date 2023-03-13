import curl
import requests
from urllib.parse import urljoin
from src.benchling import BenchlingConnection
from src.utils.base_classes import BaseConnection


class Caller:
    def __init__(self, endpoint):
        self.__setattr__('endpoint', endpoint)

    def make_request(self, method, access_token, data):
        methods = {
            "get": self.make_get,
            "post": self.make_post,
            "patch": self.make_patch,
        }

        headers = {
            'accept': 'application/json',
            'Content-Type': 'application/json',
            'Authorization': f"Bearer {access_token}"
        }

        return methods[method](headers, data)

    def make_get(self, headers, get_path):
        url = urljoin(self.__getattribute__('endpoint'), get_path)
        res = requests.get(url, headers=headers)

        if res.ok:
            print(f'Successful request. Status code: {res.status_code}.')
        else:
            print(f'Unsuccessful request. Status code: {res.status_code}. Reason: {res.reason}')
            print(f'DEBUG: {res.text}')


        return res.text

    def make_post(self, headers, json_data):
        res = requests.post(self.__getattribute__('endpoint'), json=json_data, headers=headers)
        if res.ok:
            print(f'Successful request. Status code: {res.status_code}.')
        else:
            print(f'Unsuccessful request. Status code: {res.status_code}. Reason: {res.reason}')
            print(f'DEBUG: {res.text}')

        return res

    def make_patch(self, headers, data):
        res = requests.patch(self.__getattribute__('endpoint'), json=data, headers=headers)

        if res.ok:
            print(f'Successful request. Status code: {res.status_code}.')
        else:
            print(f'Unsuccessful request. Status code: {res.status_code}. Reason: {res.reason}')
            print(f'DEBUG: {res.text}')

        return res
    
    
    def _response_handler(self, res):
        if res.ok:
            print(f'Successful request. Status code: {res.status_code}.')
        else:
            print(f'Unsuccessful request. Status code: {res.status_code}. Reason: {res.reason}')
            print(f'DEBUG: {res.text}')
            if res.status_code == '201':
                print(f'Token expired, renewing token.')
                

def export_to_service(
    json_dict: dict, 
    service_url : str,
    connection: BaseConnection,
    action : str='get', 
) -> str:

    api_caller = Caller(service_url)
    response = api_caller.make_request(action, token, json_dict)
        
    return response

def export_to_benchling(    
    json_dict: dict, 
    service_url : str,
    connection: BenchlingConnection,
    action : str='get', 
) -> str:
    if not response.ok:
        token = connection._auth_object.get_access_token()
    response = export_to_service(json_dict, service_url, connection.token, action=action)
    
    return response 
