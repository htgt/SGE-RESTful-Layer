from __future__ import annotations
import requests
from pathlib import Path
from src.utils.exceptions import NoSecretKeyException
from src.utils.base_classes import BaseClass
from time import clock_gettime, CLOCK_REALTIME
from src.rest_calls.send_calls import export_to_service
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from src.benchling import BenchlingConnection


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

        secret_path = Path(SECRET_KEY_URL)
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
        token = BenchlingToken(auth_json)
        return token
    
class BenchlingToken(BaseClass):
    def __init__(self, auth_response: dict):
        self.value = auth_response['access_token']
        self.start_time = clock_gettime(CLOCK_REALTIME)
        self.duration = auth_response['expires_in']
        self.expire_time = self.start_time + self.duration
        self.valid = True

    def check_if_expired(self):
        if self.start_time + self.duration > self.expire_time:
            return True
        else:
            return False
        

        return auth_json['access_token']

def export_to_benchling(
    json_dict: dict,
    service_url : str,
    connection: BenchlingConnection,
    action : str = 'get',
) -> str:
    # Check if token has expired.
    if connection.token.check_if_expired():
        connection.get_store_token()
    response = export_to_service(json_dict, service_url, connection.token.value, action=action)
    # If somehow it expired within timer, or token is otherwise invalid, get a new token and try again.
    response = export_to_service(json_dict, service_url, connection.token, action=action)
    if response.status_code in ["400", "401", "403"] and not response.ok:
        connection.get_store_token()
        response = export_to_service(json_dict, service_url, connection.token, action=action)

    return response
