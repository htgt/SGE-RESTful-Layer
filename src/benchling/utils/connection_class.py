from src.benchling.utils.auth_utils import APIConnector
from src.utils.base_classes import BaseClass

class BenchlingConnection(BaseClass):
    def __init__(self, client_id, benchling_tenant, secret_key, token_url):
        self.tenant = benchling_tenant
        self.secret_key = secret_key
        self.client_id = client_id
        self.token_url = token_url
        self.get_store_token()

    def get_store_token(self):
        self._auth_object = APIConnector(self.token_url, self.client_id, self.secret_key)
        if self._auth_object:
            print('BenchlingConnection initialized')
        else:
            raise (Exception("APIConnector failed to make _auth_object."))
        self.token = self._auth_object.token
    