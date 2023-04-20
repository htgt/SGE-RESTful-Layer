from src.biology.targeton_oligos import TargetonOligo
from src.benchling.utils.schemas import get_strand_dropdown_id, get_chromosome_dropdown_id
from src.core.post_targeton_oligos import send_targeton_oligo_post_request
from src.benchling import benchling_schema_ids



def post_targeton_oligos(oligo_data: dict, url: str) -> None:
    for oligo in oligo_data:
        targeton_oligo = TargetonOligo(oligo_data[oligo])
        packet = prepare_targeton_oligo_packet(oligo, targeton_oligo, benchling_schema_ids)
        response = send_targeton_oligo_post_request(packet, url)
        print(response)


def prepare_targeton_oligo_packet(name: str, targeton_oligo: TargetonOligo) -> dict:
    packet = {}
    schema_ids = benchling_schema_ids.ids

    packet['name'] = name
    packet['folderId'] = schema_ids['default_folder_id']
    packet['schemaId'] = schema_ids['schemas']['targeton_oligo_schema_id']

    packet['fields'] = as_benchling_entity(targeton_oligo)

    return packet


def as_benchling_entity(targeton_oligo: TargetonOligo) -> dict:
    entity = {
        'Action Vector' : {
            'value' : targeton_oligo.action_vector
        },
        'Ext Vector' : {
            'value' : targeton_oligo.ext_vector
        },
        'Ref. Start Position' : {
            'value' : int(targeton_oligo.ref_start)
        },
        'Ref. End Position' : {
            'value' : int(targeton_oligo.ref_end)
        },
        'R2 Start Position' : {
            'value' : int(targeton_oligo.r2_start)
        },
        'R2 End Position' : {
            'value' : int(targeton_oligo.r2_end)
        },
        'Reference Chromosome' : {
            'value' : get_chromosome_dropdown_id(targeton_oligo.ref_chromosome)
        },
        'Reference Strand' : {
            'value' : get_strand_dropdown_id(targeton_oligo.ref_strand, "plus")
        },
        'sgRNA Vector' : {
            'value' : targeton_oligo.sgrna_vector
        },
    }

    return entity
