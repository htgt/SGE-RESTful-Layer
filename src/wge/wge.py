import json
import requests

from src.domain.guideRNA import GuideRNA
from src.rest_calls.send_calls import export_to_service
from src.benchling import benchling_connection

def patch_wge_data_to_service(self, event_data : dict) -> dict:
    wge_data = query_wge_by_id(event_data['wge_id'])
    grna_class = prepare_guide_rna_class(event_data, wge_data)
    benchling_body = grna_class.as_benchling_req_body(event_data)
    patch_url = benchling_connection.sequence_url + event_data['entity_id']
    response = export_to_service(
        benchling_body,
        patch_url,
        benchling_connection.token,
        'patch',
    )
    return response

def query_wge_by_id(wge_id : str) -> dict:
    url = "https://wge.stemcell.sanger.ac.uk/api/crispr_by_id?species=Grch38&id=" + wge_id
    wge_packet = requests.get(url)

    return wge_packet.json()
     
def prepare_guide_rna_class(event_data : dict, wge_data : dict) -> GuideRNA:
    wge_id = event_data['wge_id']

    grna_data = wge_data[wge_id]
    strand = calculate_strand(grna_data['pam_right'])
    wge_link = build_wge_link(wge_id)
    species = get_wge_species(grna_data['species_id'])

    grna_dict = {
        'seq' : grna_data['seq'],
        'targeton' : event_data['targeton_id'],
        'strand' : strand,
        'wge_id' : wge_id,
        'wge_link' : wge_link,
        'off_targets' : grna_data['off_target_summary'],
        'species' : species,
    }
    print(grna_dict)
    
    grna_class = GuideRNA(grna_dict)

    return grna_class

def calculate_strand(pam : int) -> chr:
    strands = {
        0 : 'sfso_qKNl7o1M', # -
        1 : 'sfso_DqRsZ1Cg', # +
    }
    return strands[pam]

def build_wge_link(wge_id : int) -> str:
    return 'https://wge.stemcell.sanger.ac.uk/crispr/' + str(wge_id)

def get_wge_species(species_id : int) -> str:
    species = {
        1 : 'sfso_gWKuC1ge', #Grch37 - Homo Sapiens
        2 : 'sfso_gjQvG19Z', #Mouse - Mus musculus
        4 : 'sfso_gWKuC1ge', #Grch38 - Homo Sapiens
    }
    return species[species_id]

