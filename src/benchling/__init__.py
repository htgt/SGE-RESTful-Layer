from .auth_utils import APIConnector
from warnings import warn
from utils.exceptions import NoBenchlingEnvMatchWarning

CLIENT_ID = '7fd79123-bff9-4de6-9afc-81197463f016'
PROD_ENV = 'tol'  # 'prod', 'ci', 'tol'


class BenchlingConnection:
    def __init__(self):
        url = self.generate_url(prod_env=PROD_ENV)
        self.api_url = url + r'api/v2/'
        self.oligos_url = self.api_url + r'dna-oligos'
        self.blobs_url = self.api_url + r'blobs/'
        self.tasks_url = self.api_url + r'workflow-tasks/'
        self.tasks_output_url = self.api_url + r'workflow-outputs'
        self.token_url = self.api_url + r'token'

        _auth_object = APIConnector(self.token_url, CLIENT_ID)
        if _auth_object:
            print('BenchlingConnection initialized')
        else:
            raise (Exception("APIConnector failed to make _auth_object."))

        self.token = _auth_object.token

    @staticmethod
    def generate_url(prod_env='tol') -> str:
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


benchling_connection = BenchlingConnection()
