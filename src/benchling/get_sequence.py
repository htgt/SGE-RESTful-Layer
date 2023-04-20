from src.rest_calls.send_calls import Caller

import json
from urllib.parse import urljoin
from src.benchling import benchling_urls
from src.benchling.utils.request_to_benchling import request_to_benchling

def get_sequence(id):
    api_path = benchling_urls.sequence_url
    url = urljoin(api_path + "/", str(id))
    get_data = request_to_benchling(url, 'get')

    data_dict = json.loads(get_data)
    return data_dict["bases"]
