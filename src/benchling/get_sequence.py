from src.rest_calls.send_calls import Caller
from src.benchling.connection.benchling_connection import benchling_connection

import json
import posixpath
from urllib.parse import urljoin

api_path = benchling_connection.sequence_url


def get_sequence(id):
    url = urljoin(api_path + "/", str(id))
    api_caller = Caller(url)
    token = benchling_connection.token
    get_data = api_caller.make_request('get', token, url)

    data_dict = json.loads(get_data)
    return data_dict["bases"]
