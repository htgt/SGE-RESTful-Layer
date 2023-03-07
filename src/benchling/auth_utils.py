import requests
from pathlib import Path
from src.utils.exceptions import NoSecretKeyException


class APIConnector:
    def __init__(self, token_url, client_id):
        self.token_url = token_url
        self.key = self.get_secret_key()
        self.client_id = client_id
        self.auth_data = {
            "client_secret" : self.key,
            "client_id" : self.client_id,
            "grant_type" : "client_credentials"
        }
        self.token = self.get_access_token()

    def get_secret_key(self) -> str:
        # Replace with function arg and user input for url/path.
        SECRET_KEY_URL = 'src/benchling/config.cfg'
        secret_path = Path('src/benchling/config.cfg')
        if secret_path.exists():
            secret_key = open(secret_path, 'r').read()
        else:
            raise NoSecretKeyException(f"No Config file found at {SECRET_KEY_URL}")
        if len(secret_key) < 1:
            raise NoSecretKeyException(f"No secret key found at {SECRET_KEY_URL}")

        return open(SECRET_KEY_URL, 'r').read().strip('\n')

    def get_access_token(self) -> str:
        # Ideally store access token in cache with correct ttd
        # Only regenerate when cached token expires
        auth_res = requests.post(self.token_url, data=self.auth_data)
        auth_json = auth_res.json()
        
        return auth_json['access_token']
