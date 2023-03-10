from .auth_utils import APIConnector
from warnings import warn
from src.utils.exceptions import NoBenchlingEnvMatchWarning

import json
import os

# CLIENT_ID = '7fd79123-bff9-4de6-9afc-81197463f016' # tol
CLIENT_ID = '7df4bb27-81bc-4be8-b08c-afac5609a195'
PROD_ENV = 'ci'  # 'prod', 'ci', 'tol'

BENCHLING_IDS_URL = 'benchling_schema_ids.json'


class BenchlingConnection:
    def __init__(self):
        url = self.generate_url(prod_env=PROD_ENV)
        self.api_url = url + r'api/v2/'
        self.blobs_url = self.api_url + r'blobs/'
        self.oligos_url = self.api_url + r'dna-oligos'
        self.sequence_url = self.api_url + r'dna-sequence'
        self.tasks_url = self.api_url + r'workflow-tasks/'
        self.tasks_output_url = self.api_url + r'workflow-outputs'
        self.custom_entity_url = self.api_url + r'custom-entities'

        self.token_url = self.api_url + r'token'

        _auth_object = APIConnector(self.token_url, CLIENT_ID)
        if _auth_object:
            print('BenchlingConnection initialized')
        else:
            raise (Exception("APIConnector failed to make _auth_object."))

        self.token = _auth_object.token

    @staticmethod
    def generate_url(prod_env='ci') -> str:
        url = r'https://'
        prod_env_dict = {
            "tol" : r"tol-sangertest.",
            "prod" : r"ci-sanger.",
            "ci" : r"ci-sanger-test."
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
        self.ids = json.load(open(BENCHLING_IDS_URL))


benchling_connection = BenchlingConnection()
benchling_schema_ids = BenchlingSchemaIds()
