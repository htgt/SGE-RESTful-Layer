from .auth_utils import APIConnector
from warnings import warn
from src.utils.exceptions import NoBenchlingEnvMatchWarning
from src import BENCHLING_SECRET_KEY, BENCHLING_TENANT

import json


# PROD_ENV = 'ci'  # 'prod', 'ci', 'tol'
# TEST (test)
CI_TEST_CLIENT_ID = '7df4bb27-81bc-4be8-b08c-afac5609a195'
CI_TEST_BENCHLING_IDS_URL = r'schemas/ci_sanger_test_ids.json'
# MAVE (prod)
MAVE_SANGER_CLIENT_ID = 'TEMP'
MAVE_SANGER_BENCHLING_IDS_URL = r'schemas/mave_sanger_ids.json'


if BENCHLING_TENANT == 'prod':
    benchling_ids_url = MAVE_SANGER_BENCHLING_IDS_URL
    client_id = MAVE_SANGER_CLIENT_ID
elif BENCHLING_TENANT == 'test':
    benchling_ids_url = CI_TEST_BENCHLING_IDS_URL
    client_id = CI_TEST_CLIENT_ID


class BenchlingConnection:
    def __init__(self):
        self.prod_env = BENCHLING_TENANT
        self.secret_key = BENCHLING_SECRET_KEY
        url = self.generate_url(prod_env=BENCHLING_TENANT)
        self.api_url = url + r'api/v2/'
        self.blobs_url = self.api_url + r'blobs/'
        self.oligos_url = self.api_url + r'dna-oligos'
        self.sequence_url = self.api_url + r'dna-sequences'
        self.tasks_url = self.api_url + r'workflow-tasks/'
        self.tasks_output_url = self.api_url + r'workflow-outputs'
        self.custom_entity_url = self.api_url + r'custom-entities'

        self.token_url = self.api_url + r'token'

        self.get_store_token()

    def get_store_token(self):
        self._auth_object = APIConnector(self.token_url, client_id, self.secret_key)
        if self._auth_object:
            print('BenchlingConnection initialized')
        else:
            raise (Exception("APIConnector failed to make _auth_object."))
        self.token = self._auth_object.token

    @staticmethod
    def generate_url(prod_env='ci') -> str:
        url = r'https://'
        prod_env_dict = {
            "tol" : r"tol-sangertest.",
            "prod" : r"ci-sanger.",
            "test" : r"ci-sanger-test."
        }
        if prod_env in prod_env_dict:
            url = url + prod_env_dict[prod_env]
        else:
            warn(
                f"Selected benchling environment doesn't match {prod_env_dict.keys()}\nUsing {prod_env_dict['tol']}",
                NoBenchlingEnvMatchWarning)
            url = url + prod_env_dict["tol"]

        return url + r"benchling.com/"


class BenchlingSchemaIds:
    def __init__(self):
        self.ids = json.load(open(benchling_ids_url))


benchling_connection = BenchlingConnection()
benchling_schema_ids = BenchlingSchemaIds()
