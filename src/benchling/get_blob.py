from src.rest_calls.send_calls import Caller

import json
import posixpath
from urllib.parse import urljoin




def get_blob_url(id, api_path, token):
    path = posixpath.join(id, 'download-url')
    url = urljoin(api_path, path)

    api_caller = Caller(url)

    get_data = api_caller.make_request('get', token, url)

    return json.loads(get_data)["downloadURL"]
