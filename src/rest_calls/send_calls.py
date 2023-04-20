import requests
from src.utils.exceptions import ConvertToJsonError
import curl
from typing import Tuple
from urllib.parse import urljoin


class Caller:
    def __init__(self, endpoint):
        self.endpoint = endpoint

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

    def make_get(self, headers, data):
        url = urljoin(self.endpoint, data)
        res = requests.get(url, headers=headers)
        self._response_handler(res)

        return res.text

    def make_post(self, headers, data):
        res = requests.post(self.endpoint, json=data, headers=headers)
        self._response_handler(res)

        return res

    def make_patch(self, headers, data):
        res = requests.patch(self.endpoint, json=data, headers=headers)
        self._response_handler(res)

        return res

    def _response_handler(self, res):
        if res.ok:
            print(f'Successful request. Status code: {res.status_code}.')
        else:
            print(f'Unsuccessful request. Status code: {res.status_code}. Reason: {res.reason}')
            print(f'DEBUG: {res.text}')
            print(curl.parse(res))


def request_to_service(
    service_url: str,
    token: str,
    action: str,
    data: Tuple[dict, str]
) -> requests.Response:

    api_caller = Caller(service_url)
    response = api_caller.make_request(action, token, data)

    return response


def check_response_object(response_object) -> dict:
    try:
        json_response = response_object.json()

    except Exception as err:
        print(err)
        raise ConvertToJsonError(err)

    return json_response
