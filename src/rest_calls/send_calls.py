import requests
from benchling.auth_uitls import AuthUtils


class Caller:
    def __init__(self, endpoint):
        self.__setattr__('endpoint', endpoint)

    def make_post(self, json_data):

        access_token = AuthUtils.get_access_token()

        print(access_token)

        headers = {'Authorization': f"Bearer {access_token}"}

        res = requests.post(self.__getattribute__('endpoint'), json=json_data, headers=headers)

        if res.ok:
            print(f'Successful request. Status code: {res.status_code}.')
        else:
            print(f'Unsuccessful request. Status code: {res.status_code}. Reason: {res.reason}')
            print(f'DEBUG: {res.text}')

        return res

