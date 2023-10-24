from src.rest_calls.send_calls import Caller

import json
from urllib.parse import urljoin
from src.benchling import benchling_urls
from src.benchling.utils.request_to_benchling import request_to_benchling

def get_sequence(id):
    api_path = benchling_urls.sequence_url
    data_dict = request_to_benchling(api_path, 'get', str(id)).json()

    return data_dict["bases"]
