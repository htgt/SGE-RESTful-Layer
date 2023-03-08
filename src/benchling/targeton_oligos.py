import json

from src.domain.targeton_oligos import TargetonOligo
from src.rest_calls.send_calls import export_to_service
# ts_UP3gzN10

def post_targeton_oligos(oligo_data: dict):
    for oligo in oligo_data:
        targeton_oligo = TargetonOligo(oligo_data[oligo])
        #print(targeton_oligo.as_benchling_entity())
        prepare_targeton_oligo_packet(targeton_oligo)

def prepare_targeton_oligo_packet(targeton_oligo: TargetonOligo):
    packet = {}
    packet['fields'] = as_benchling_entity(targeton_oligo)

    return packet


def as_benchling_entity(targeton_oligo: TargetonOligo) -> dict:
    entity = {
        'Action Vector' : targeton_oligo.action_vector,
        'Ext Vector' : targeton_oligo.ext_vector,
        'Ref. Start Position' : targeton_oligo.ref_start,
        'Ref. End Position' : targeton_oligo.ref_end,
        'Ref2. Start Position' : targeton_oligo.r2_start,
        'Ref2. End Position' : targeton_oligo.r2_end,
        'Reference Chromosome' : targeton_oligo.ref_chromosome,
        'Reference Strand' : targeton_oligo.ref_strand,
        'Version' : targeton_oligo.version,
        'sgRNA Vector' : targeton_oligo.sgrna_vector,
    }

    return entity
