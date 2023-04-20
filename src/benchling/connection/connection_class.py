from src.benchling.connection.auth_utils import APIConnector
from warnings import warn
from src.utils.exceptions import NoBenchlingEnvMatchWarning
from src.utils.base_classes import BaseClass

class BenchlingConnection(BaseClass):
    def __init__(self, client_id, benchling_tenant, secret_key):
        self.tenant = benchling_tenant
        self.secret_key = secret_key
        self.client_id = client_id
        self.get_store_token()

    def get_store_token(self):
        self._auth_object = APIConnector(self.token_url, self.client_id, self.secret_key)
        if self._auth_object:
            print('BenchlingConnection initialized')
        else:
            raise (Exception("APIConnector failed to make _auth_object."))
        self.token = self._auth_object.token
    
class BenchlingUrls(BaseClass):
    def __init__(self, tenant) -> None:
        url = self.generate_url(tenant)
        self.api_url = url + r'api/v2/'
        self.blobs_url = self.api_url + r'blobs/'
        self.oligos_url = self.api_url + r'dna-oligos'
        self.sequence_url = self.api_url + r'dna-sequences'
        self.tasks_url = self.api_url + r'workflow-tasks/'
        self.tasks_output_url = self.api_url + r'workflow-outputs'
        self.custom_entity_url = self.api_url + r'custom-entities'
        self.token_url = self.api_url + r'token'
        
    @staticmethod
    def generate_url(tenant: str ='ci') -> str:
        url = r'https://'
        tenant_dict = {
            "tol" : r"tol-sangertest.",
            "prod" : r"mave-sanger.",
            "test" : r"ci-sanger-test.",
            "unittest" : r"unittest."
        }
        if tenant in tenant_dict:
            url = url + tenant_dict[tenant]
        else:
            warn(
                f"Selected benchling environment {tenant} doesn't match {tenant_dict.keys()}\nUsing {tenant_dict['tol']}",
                NoBenchlingEnvMatchWarning)
            url = url + tenant_dict["tol"]
        
        url = url + r"benchling.com/"
        
        return url
    