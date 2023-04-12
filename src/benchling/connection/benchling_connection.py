from src import BENCHLING_SECRET_KEY, BENCHLING_TENANT
from src.benchling import client_id
from src.benchling.connection.connection_class import BenchlingConnection

benchling_connection = BenchlingConnection(client_id, BENCHLING_TENANT, BENCHLING_SECRET_KEY)