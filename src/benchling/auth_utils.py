from __future__ import annotations
import requests


class APIConnector:
    def __init__(self, token_url, client_id, secret_key):
        self.token_url = token_url
        self.key = secret_key
        self.client_id = client_id
        self.auth_data = {
            "client_secret" : self.key,
            "client_id" : self.client_id,
            "grant_type" : "client_credentials"
        }
        self.token = self.get_access_token()

    def get_access_token(self) -> str:
        # Ideally store access token in cache with correct ttd
        # Only regenerate when cached token expires
        print('APIConnector get new token')

        auth_res = requests.post(self.token_url, data=self.auth_data)
        auth_json = auth_res.json()

        result = auth_json
        if 'access_token' in auth_json:
            result = auth_json['access_token']

        return result
