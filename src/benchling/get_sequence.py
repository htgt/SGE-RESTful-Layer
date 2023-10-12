from src.rest_calls.send_calls import Caller

import json
from urllib.parse import urljoin
from src.benchling import benchling_urls
from src.benchling.utils.request_to_benchling import request_to_benchling

def get_sequence(id):
    print('Get sequence from Benchling')

    api_path = benchling_urls.sequence_url
    get_data = request_to_benchling(api_path, 'get', str(id))

    print('GET SEQUENCE RESULT::::', get_data)
    data_dict = json.loads(get_data)
    print('DATA DICT::::', data_dict)

    return data_dict["bases"]
