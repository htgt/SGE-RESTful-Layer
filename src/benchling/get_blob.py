from src.rest_calls.send_calls import Caller
from src.benchling import benchling_connection

import json
import posixpath
from urllib.parse import urljoin

api_path = 'https://tol-sangertest.benchling.com/api/v2/blobs/'


def get_blob_url(id):
    path = posixpath.join(id, 'download-url')

    url = urljoin(api_path, path)

    api_caller = Caller(url)
    token = benchling_connection.token

    get_data = api_caller.make_request('get', token, url)

    return json.loads(get_data)["downloadURL"]
