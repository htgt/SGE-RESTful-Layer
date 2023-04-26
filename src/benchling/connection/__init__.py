from src import BENCHLING_SECRET_KEY, BENCHLING_TENANT
from src.benchling import client_id, benchling_urls
from src.benchling.utils.connection_class import BenchlingConnection

print("Initialising connection...")
benchling_connection = BenchlingConnection(client_id, BENCHLING_TENANT, BENCHLING_SECRET_KEY, benchling_urls.token_url)