import curl
import requests
from urllib.parse import urljoin
from src.benchling import BenchlingConnection


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
        # print(url)
        res = requests.get(url, headers=headers)

        if res.ok:
            print(f'Successful request. Status code: {res.status_code}.')
        else:
            print(f'Unsuccessful request. Status code: {res.status_code}. Reason: {res.reason}')
            print(f'DEBUG: {res.text}')

        # print(res)

        return res.text
    # res.json()

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
            print(curl.parse(res))

        return res

def export_to_benchling(json_dict: dict, benchling_connection: BenchlingConnection) -> str:
    api_caller = Caller(benchling_connection.oligos_url)
    token = benchling_connection.token

    try:
        oligos_id = api_caller.make_request('post', token, json_dict).json()['id']

    except Exception as err:
        raise Exception(err)

    return oligos_id
