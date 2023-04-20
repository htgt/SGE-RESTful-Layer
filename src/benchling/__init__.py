from typing import Tuple
from src import BENCHLING_TENANT
import json
from warnings import warn
from src.utils.exceptions import NoBenchlingEnvMatchWarning
from src.benchling.connection.connection_class import BenchlingUrls

# TENANT = 'ci'  # 'prod', 'ci', 'tol', 'test', 'unittest'
# TEST (test)
CI_TEST_CLIENT_ID = '7df4bb27-81bc-4be8-b08c-afac5609a195'
CI_TEST_BENCHLING_IDS_URL = r'schemas/ci_sanger_test_ids.json'
# MAVE (prod)
MAVE_SANGER_CLIENT_ID = 'a669776b-16f2-431c-944a-3e01318a14a6'
MAVE_SANGER_BENCHLING_IDS_URL = r'schemas/mave_sanger_ids.json'
# UNITTEST (unittest)
UNITTEST_CLIENT_ID = 'unittest'
UNITTEST_BENCHLING_IDS_URL = r'schemas/ci_sanger_test_ids.json'

class BenchlingSchemaIds:
    def __init__(self, benchling_ids_url):
        with open(benchling_ids_url) as f:
            self.ids = json.load(f)

    
def get_tenant_ids(tenant:str) -> Tuple[str, str]:
    if tenant == 'prod':
        benchling_ids_url = MAVE_SANGER_BENCHLING_IDS_URL
        client_id = MAVE_SANGER_CLIENT_ID
    elif tenant == 'test':
        benchling_ids_url = CI_TEST_BENCHLING_IDS_URL
        client_id = CI_TEST_CLIENT_ID
    else:
        warn(f"No matching benchling tenant found {tenant}, using test.", NoBenchlingEnvMatchWarning)
        benchling_ids_url = CI_TEST_BENCHLING_IDS_URL
        client_id = CI_TEST_CLIENT_ID
        
    return client_id, benchling_ids_url

client_id, benchling_ids_url = get_tenant_ids(BENCHLING_TENANT)
benchling_schema_ids = BenchlingSchemaIds(benchling_ids_url)
benchling_urls = BenchlingUrls(BENCHLING_TENANT)
