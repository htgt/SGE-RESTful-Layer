import requests
from urllib.parse import urljoin

class Caller:
    def __init__(self, endpoint):
        self.__setattr__('endpoint', endpoint)

    def make_request(self, method, access_token, data):
        methods = {
            "get": self.make_get,
            "post": self.make_post,
            "patch": self.make_patch
        }

        headers = {'Authorization': f"Bearer {access_token}"}

        return methods[method]( headers, data)

    def make_get(self, headers, get_path):
        url = urljoin(self.__getattribute__('endpoint'), get_path)

        res = requests.get(url, headers=headers)

        if res.ok:
            print(f'Successful request. Status code: {res.status_code}.')
        else:
            print(f'Unsuccessful request. Status code: {res.status_code}. Reason: {res.reason}')
            print(f'DEBUG: {res.text}')

        print(res)

        return res.text

    def make_post(self, headers, json_data):
        res = requests.post(self.__getattribute__('endpoint'), json=json_data, headers=headers)

        if res.ok:
            print(f'Successful request. Status code: {res.status_code}.')
        else:
            print(f'Unsuccessful request. Status code: {res.status_code}. Reason: {res.reason}')
            print(f'DEBUG: {res.text}')

        return res

    def make_patch(self, headers, params_data):
        res = requests.patch(self.__getattribute__('endpoint'), params=params_data, headers=headers)

        if res.ok:
            print(f'Successful request. Status code: {res.status_code}.')
        else:
            print(f'Unsuccessful request. Status code: {res.status_code}. Reason: {res.reason}')
            print(f'DEBUG: {res.text}')

        return res
