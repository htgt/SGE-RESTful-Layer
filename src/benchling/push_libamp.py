from src.benchling.connection.benchling_connection import BenchlingConnection, benchling_connection
from src.benchling import benchling_schema_ids
from src.benchling.utils.export_to_benchling import export_to_benchling
from src.benchling.archive_entity import archive_oligo
from src.biology.libamp_primers import LibampPrimer



def primer_to_benchling_json(primer: LibampPrimer, ids) -> dict:
    return {
        "bases": primer.sequence,
        "name": primer.name,
        "fields": {
            "GC Content (%)": {
                "value": primer.gc_content,
            },
            "Genome Location": {
                "value": primer.chr_start,
            },
            "Primer Direction (short form)": {
                "value": ids["dropdowns"]["forward_direction"] if primer.strand == "left" else ids["dropdowns"]["reverse_direction"],
            },
            "LibAmp Primer Type": {
                "value":  ids["dropdowns"]["libamp_forward"] if primer.strand == "left" else ids["dropdowns"]["libamp_reverse"],
            },
            "Primer Score": {
                "value": primer.score,
            },
            "Product Size (bp)":{
                "value": primer.product_size,
            },
            "Tm (°C)": {
                "value": primer.melting_temp,
            },
        },
        "folderId": ids["default_folder_id"],
        "schemaId": ids["schemas"]["libamp_schema_id"]
    }


def call_export_primer_pair(primer_left: LibampPrimer, primer_right: LibampPrimer):
    url = benchling_connection.oligos_url

    return export_primer_pair(
        primer_left,
        primer_right,
        export_to_benchling,
        archive_oligo,
        url,
        benchling_connection
    )

def export_primer_pair(
        primer_left: LibampPrimer,
        primer_right: LibampPrimer,
        export_function,
        archive_function,
        url: str,
        connection: BenchlingConnection,
) -> list:

    primer_left_json = primer_to_benchling_json(primer_left, benchling_schema_ids.ids)
    primer_right_json = primer_to_benchling_json(primer_right, benchling_schema_ids.ids)

    left_response = export_function(
        primer_left_json,
        url,
        connection,
        'post',
    )

    if left_response.ok:
        left_primer_id =  left_response.json()["id"]

        right_response = export_function(
            primer_right_json,
            url,
            connection,
            'post',
        )

        if not right_response.ok:
            archive_function(left_primer_id)

            raise Exception(right_response)

    else:
        raise Exception(left_response)

    return [left_response.json(), right_response.json()]
