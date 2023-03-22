from __future__ import annotations
import requests
from src.utils.exceptions import NoSecretKeyException
from src.rest_calls.send_calls import export_to_service, check_response_object
from typing import TYPE_CHECKING
from dotenv import load_dotenv
import os

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
        secret_key = os.getenv('BENCHLING_SECRET_KEY')
        if not secret_key:
            try:
                load_dotenv(".env")
                secret_key = os.getenv('BENCHLING_SECRET_KEY')
                if not secret_key:
                    raise NoSecretKeyException(f"Environmental variable Benchling secret key is empty.")
            except:
                raise NoSecretKeyException(f"No secret key found in Enviromental variables or .env")

        if len(secret_key) < 1:
            raise NoSecretKeyException(f"Empty secret key stored in .env")
        return secret_key

    def get_access_token(self) -> str:
        # Ideally store access token in cache with correct ttd
        # Only regenerate when cached token expires
        auth_res = requests.post(self.token_url, data=self.auth_data)
        auth_json = auth_res.json()

        return auth_json['access_token']

def export_to_benchling(
    json_dict: dict,
    service_url : str,
    connection: BenchlingConnection,
    action : str = 'get',
) -> str:

    response = export_to_service(json_dict, service_url, connection.token, action=action)
    if response.status_code in ["400", "401", "403"] and not response.ok:
        connection.get_store_token()
        response = export_to_service(json_dict, service_url, connection.token, action=action)
    
    json_response = check_response_object(response)

    return json_response
