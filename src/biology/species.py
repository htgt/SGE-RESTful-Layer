def get_species_name_by_id(id):
    return SPECIES[id]["name"]

SPECIES = {
    'Grch37': {
        'name': 'homo_sapiens',
    },
    'Grch38': {
        'name': 'homo_sapiens',
    },
    'Mouse': {
        'name': 'mus_musculus',
    },
}