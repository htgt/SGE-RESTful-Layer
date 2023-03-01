import json
import requests

from src.domain.guideRNA import GuideRNA

def query_wge_by_id(wge_id : str) -> dict:
    url = "https://wge.stemcell.sanger.ac.uk/api/crispr_by_id?species=Grch38&id=" + wge_id
    wge_packet = requests.get(url)

    return wge_packet.json()
     
def prepare_guide_rna_entity(event_data : dict, wge_data : dict) -> dict:
    #{'1107691594': {'seq': 'AGGATTTTGGTAATTGGGTTTGG', 'genic': 1, 'chr_name': '13', 'off_target_summary': '{0: 1, 1: 0, 2: 1, 3: 15, 4: 204}', 'chr_start': 42958750, 'chr_end': 42958772, 'pam_right': 1, 'exonic': 0, 'species_id': 4, 'id': 1107691594}}
    wge_id = event_data['wge_id']

    grna_data = wge_data[wge_id]
    strand = calculate_strand(grna_data['pam_right'])
    wge_link = build_wge_link(wge_id)
    species = get_wge_species(grna_data['species_id'])

    grna_entity = {
        'Guide Sequence' : grna_data['seq'],
        'Targeton' : event_data['targeton_id'],
        'Strand' : strand,
        'WGE ID' : wge_id,
        'WGE Hyperlink' : wge_link,
        'Off Target Summary Data' : grna_data['off_target_summary'],
        'Species' : species,
    }

    print(grna_entity)
    return grna_entity

def post_guide_rna_to_benchling(event_data, grna_entity):
    return

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
