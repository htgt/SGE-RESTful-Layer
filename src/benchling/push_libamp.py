from . import benchling_connection, benchling_schema_ids
from src.rest_calls.send_calls import export_to_service


def primer_to_benchling_json(primer, ids) -> dict:
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

def export_primer_pair(primer_left, primer_right) -> list:
    url = benchling_connection.oligos_url
    token = benchling_connection.token

    primer_left_json = primer_to_benchling_json(primer_left, benchling_schema_ids.ids)
    primer_right_json = primer_to_benchling_json(primer_right, benchling_schema_ids.ids)

    try:
        left_response = export_to_service(
            primer_left_json,
            url,
            token,
            'post',
        )
        right_response= export_to_service(
            primer_right_json,
            url,
            token,
            'post',
        )

        print("Result: ", left_response, right_response)

    except Exception as err:
        raise Exception(err)

    return [left_response, right_response]
