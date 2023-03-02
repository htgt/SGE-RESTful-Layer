import json
import requests

from src.domain.guideRNA import GuideRNA

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
        0 : '-',
        1 : '+',
    }
    return strands[pam]

def build_wge_link(wge_id : int) -> str:
    return 'https://wge.stemcell.sanger.ac.uk/crispr/' + str(wge_id)

def get_wge_species(species_id : int) -> str:
    species = {
        1 : 'Human',
        2 : 'Mouse',
        3 : 'Pig',
        4 : 'Grch38'
    }
    return species[species_id]
