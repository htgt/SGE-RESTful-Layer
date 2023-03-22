import requests

from src.biology.guideRNA import GuideRNA
from src.rest_calls.send_calls import export_to_service_json_response
from src.benchling import benchling_connection
from src.benchling.utils.schemas import get_strand_dropdown_id


def query_wge_by_id(wge_id : str) -> dict:
    url = "https://wge.stemcell.sanger.ac.uk/api/crispr_by_id?species=Grch38&id=" + str(wge_id)
    wge_packet = requests.get(url)

    return wge_packet.json()


def prepare_guide_rna_class(event_data : dict, wge_data : dict) -> GuideRNA:
    wge_id = event_data['wge_id']

    grna_data = wge_data[wge_id]
    #strand = get_strand_dropdown_id(grna_data['pam_right'])
    wge_link = build_wge_link(wge_id)
    species = get_wge_species(grna_data['species_id'])

    grna_dict = {
        'seq' : grna_data['seq'],
        'targeton' : event_data['targeton_id'],
        #'strand' : strand,
        'wge_id' : wge_id,
        'wge_link' : wge_link,
        'off_targets' : grna_data['off_target_summary'],
        'species' : species,
    }

    grna_class = GuideRNA(grna_dict)

    return grna_class


def build_wge_link(wge_id : int) -> str:
    return 'https://wge.stemcell.sanger.ac.uk/crispr/' + str(wge_id)


def get_wge_species(species_id : int) -> str:
    species = {
        1 : 'Grch37', # Homo Sapiens
        2 : 'Mouse', # Mus musculus
        4 : 'Grch38', # Homo Sapiens
    }
    return species[species_id]


def transform_wge_event(data):
    data_entity = data['detail']['entity']

    wge_grna_data = {}
    wge_grna_data['folder_id'] = data_entity['folderId']
    wge_grna_data['entity_id'] = data_entity['id']
    wge_grna_data['wge_id'] = data_entity['fields']['WGE ID']['value']
    wge_grna_data['targeton_id'] = data_entity['fields']['Targeton']['value']
    wge_grna_data['schema_id'] = data_entity['schema']['id']
    wge_grna_data['name'] = data_entity['name']

    return wge_grna_data
