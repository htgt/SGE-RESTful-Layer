import json

from src.domain.targeton_oligos import TargetonOligo
from src.rest_calls.send_calls import export_to_service
from src.benchling import benchling_schema_ids, benchling_connection
from src.utils.schemas import get_strand_dropdown_id, get_chromosome_dropdown_id

def post_targeton_oligos(oligo_data: dict):
    for oligo in oligo_data:
        targeton_oligo = TargetonOligo(oligo_data[oligo])
        packet = prepare_targeton_oligo_packet(oligo, targeton_oligo)
        send_targeton_oligo_post_request(packet)


def prepare_targeton_oligo_packet(name: str, targeton_oligo: TargetonOligo):
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
        'Ref2. Start Position' : {
            'value' : int(targeton_oligo.r2_start)
        },
        'Ref2. End Position' : {
            'value' : int(targeton_oligo.r2_end)
        },
        'Reference Chromosome' : {
            'value' : get_chromosome_dropdown_id(targeton_oligo.ref_chromosome)
        },
        'Reference Strand' : {
            'value' : get_strand_dropdown_id(targeton_oligo.ref_strand)
        },
        'sgRNA Vector' : {
            'value' : targeton_oligo.sgrna_vector
        },
    }

    return entity


def send_targeton_oligo_post_request(body: dict):
    response = export_to_service(
        body,
        benchling_connection.custom_entity_url,
        benchling_connection.token,
        'post'
    )
    print(response)
