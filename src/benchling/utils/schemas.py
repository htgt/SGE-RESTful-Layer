from src.benchling import benchling_schema_ids


def get_strand_dropdown_id(sign : chr, type: str = "sense") -> str:
    dropdowns = benchling_schema_ids.ids['dropdowns']

    if type == "sense":
        fwd = dropdowns['sense']
        rev = dropdowns['antisense']
    elif type == "plus":
        fwd = dropdowns['plus']
        rev = dropdowns['minus']

    strands = {
        0 : rev,
        1 : fwd,
        '+' : fwd,
        '-' : rev,
    }

    return strands[sign]

def get_chromosome_dropdown_id(chr_id: str) -> str:
    chromosomes = benchling_schema_ids.ids['dropdowns']['chromosomes']
    return chromosomes[chr_id]
