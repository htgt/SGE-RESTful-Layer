import requests

class AuthUtils:
    def __init__(self):
        self.token_url = 'https://tol-sangertest.benchling.com/api/v2/token'
        self.key = self.get_secret_key()
        self.client_id = "7fd79123-bff9-4de6-9afc-81197463f016"
        self.auth_data = {"client_secret": self.key,
            "client_id"             : self.client_id,
            "grant_type"            : "client_credentials"}
        self.token = self.get_access_token()

    def get_secret_key(self):
        SECRET_KEY_URL = 'src/benchling/config.cfg'

        return open(SECRET_KEY_URL, 'r').read()

    def get_access_token(self):
        # Ideally store access token in cache with correct ttd
        # Only regenerate when cached token expires
        auth_res = requests.post(self.token_url, data=self.auth_data)
        auth_json = auth_res.json()

        return auth_json['access_token']
