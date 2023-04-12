from src.rest_calls.send_calls import Caller

import json
from urllib.parse import urljoin

def get_sequence(id, api_path, token):
    url = urljoin(api_path + "/", str(id))
    api_caller = Caller(url)
    get_data = api_caller.make_request('get', token, url)

    data_dict = json.loads(get_data)
    return data_dict["bases"]
