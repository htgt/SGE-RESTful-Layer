from src.benchling.utils.request_to_benchling import request_to_benchling
from src.benchling import benchling_urls
import json
import posixpath
from urllib.parse import urljoin




def get_blob_url(id):
    api_path = benchling_urls.blobs_url
    path = posixpath.join(id, 'download-url')
    url = urljoin(api_path, path)

    get_data = request_to_benchling(url, 'get')

    return json.loads(get_data)["downloadURL"]
