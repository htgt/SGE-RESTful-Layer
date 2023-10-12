from src.rest_calls.send_calls import Caller

import json
from urllib.parse import urljoin
from src.benchling import benchling_urls
from src.benchling.utils.request_to_benchling import request_to_benchling

def get_sequence(id):
    print('Get sequence from Benchling')

    path = urljoin(benchling_urls.sequence_url, str(id))
    print('Sequence URL:::::', path)
    get_data = request_to_benchling(path, 'get')

    data_dict = json.loads(get_data)
    return data_dict["bases"]
