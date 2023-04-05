from .auth_utils import APIConnector
from warnings import warn
from src.utils.exceptions import NoBenchlingEnvMatchWarning
from src import BENCHLING_SECRET_KEY, BENCHLING_TENANT
from typing import Tuple

import json

# TENANT = 'ci'  # 'prod', 'ci', 'tol'
# TEST (test)
CI_TEST_CLIENT_ID = '7df4bb27-81bc-4be8-b08c-afac5609a195'
CI_TEST_BENCHLING_IDS_URL = r'schemas/ci_sanger_test_ids.json'
# MAVE (prod)
MAVE_SANGER_CLIENT_ID = 'a669776b-16f2-431c-944a-3e01318a14a6'
MAVE_SANGER_BENCHLING_IDS_URL = r'schemas/mave_sanger_ids.json'

class BenchlingConnection:
    def __init__(self, client_id, benchling_tenant):
        self.tenant = benchling_tenant
        self.secret_key = BENCHLING_SECRET_KEY
        url = self.generate_url(tenant=benchling_tenant)
        self.api_url = url + r'api/v2/'
        self.blobs_url = self.api_url + r'blobs/'
        self.oligos_url = self.api_url + r'dna-oligos'
        self.sequence_url = self.api_url + r'dna-sequences'
        self.tasks_url = self.api_url + r'workflow-tasks/'
        self.tasks_output_url = self.api_url + r'workflow-outputs'
        self.custom_entity_url = self.api_url + r'custom-entities'
        self.token_url = self.api_url + r'token'
        self.client_id = client_id
        self.get_store_token()

    def get_store_token(self):
        self._auth_object = APIConnector(self.token_url, self.client_id, self.secret_key)
        if self._auth_object:
            print('BenchlingConnection initialized')
        else:
            raise (Exception("APIConnector failed to make _auth_object."))
        self.token = self._auth_object.token

    @staticmethod
    def generate_url(tenant='ci') -> str:
        url = r'https://'
        tenant_dict = {
            "tol" : r"tol-sangertest.",
            "prod" : r"mave-sanger.",
            "test" : r"ci-sanger-test."
        }
        if tenant in tenant_dict:
            url = url + tenant_dict[tenant]
        else:
            warn(
                f"Selected benchling environment doesn't match {tenant_dict.keys()}\nUsing {tenant_dict['tol']}",
                NoBenchlingEnvMatchWarning)
            url = url + tenant_dict["tol"]

        return url + r"benchling.com/"


class BenchlingSchemaIds:
    def __init__(self, benchling_ids_url):
        self.ids = json.load(open(benchling_ids_url))

    
def get_tenant_ids(tenant:str) -> Tuple[str, str]:
    if tenant == 'prod':
        benchling_ids_url = MAVE_SANGER_BENCHLING_IDS_URL
        client_id = MAVE_SANGER_CLIENT_ID
    elif tenant == 'test':
        benchling_ids_url = CI_TEST_BENCHLING_IDS_URL
        client_id = CI_TEST_CLIENT_ID
        
    return client_id, benchling_ids_url


client_id, benchling_ids_url = get_tenant_ids(BENCHLING_TENANT)
benchling_schema_ids = BenchlingSchemaIds(benchling_ids_url)
benchling_connection = BenchlingConnection(client_id, BENCHLING_TENANT)

