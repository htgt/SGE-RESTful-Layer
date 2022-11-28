from .auth_utils import APIConnector

class BenchlingConnection:
    def __init__(self):
        _auth_object = APIConnector()

        print('BenchlingConnection initialized')
        self.token = _auth_object.token

benchling_connection = BenchlingConnection()