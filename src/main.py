from src.benchling import create_oligos, create_gRNA

data = {
    'id': '1059134679',
    'sequence': 'CAGTAGACACATGGTATTG',
    'gene_name': 'A1CF'
}

grna_from_benchling = create_gRNA.export_grna_to_benchling(data)

guide_rna_event_data = grna_from_benchling
create_oligos.export_oligos_to_benchling(guide_rna_event_data)
