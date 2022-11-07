from src.benchling.auth_uitls import AuthUtils
from src.rest_calls.send_calls import Caller

import posixpath
from urllib.parse import urljoin

api_path = 'https://tol-sangertest.benchling.com/api/v2/blobs/'

def get_blob_url(id):
    path = posixpath.join(id, 'download-url')

    url = urljoin(api_path, path)

    api_caller = Caller(url)
    auth_object = AuthUtils()
    token = auth_object.token

    get_data = api_caller.make_request('get', token, url)

    return get_data