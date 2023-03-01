import json
import requests

def query_wge_by_id(wge_id : str) -> dict:
    url = "https://wge.stemcell.sanger.ac.uk/api/crispr_by_id?species=Grch38&id=" + wge_id
    wge_packet = requests.get(url)
    print(wge_packet)
    return wge_packet.json()
