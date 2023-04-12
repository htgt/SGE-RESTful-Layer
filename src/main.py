from src.benchling import create_oligos, create_gRNA, benchling_connection

data = {
    'id': '1059134679',
    'sequence': 'CAGTAGACACATGGTATTG',
    'gene_name': 'A1CF'
}

create_gRNA.export_grna_to_benchling(data, benchling_connection)

guide_rna_event_data = grna_from_benchling
create_oligos.export_oligos_to_benchling(guide_rna_event_data)
