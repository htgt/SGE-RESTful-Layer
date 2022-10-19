import requests


class Caller:

    def get_access_token(self):
        secret_key = open('config.cfg', 'r').read()

        # Ideally store access token in cache with correct ttd
        # Only regenerate when cached token expires
        auth_data = {
            "client_secret": secret_key,
            "client_id": "a5f3e7b7-fef4-4c00-8249-f75c39c686f7",
            "grant_type": "client_credentials"
        }
        auth_res = requests.post('https://tol-sangertest.benchling.com/api/v2/token', data=auth_data)
        auth_json = auth_res.json()

        return auth_json['access_token']

    def make_post(self, endpoint, json_data):

        access_token = self.get_access_token()

        print(access_token)

        headers = {'Authorization': f"Bearer {access_token}"}

        res = requests.post(endpoint, json=json_data, headers=headers)

        if res.ok:
            print(f'Successful request. Status code: {res.status_code}.')
        else:
            print(f'Unsuccessful request. Status code: {res.status_code}. Reason: {res.reason}')
            print(f'DEBUG: {res.text}')

        return res

