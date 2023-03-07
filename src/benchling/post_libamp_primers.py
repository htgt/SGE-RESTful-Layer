from src.rest_calls.send_calls import Caller
from . import benchling_connection

from src.domain.libamp_primers import LibampPrimer
import json
import sys
sys.path.append("..")


def post_libamp_primers(data):
    for primer_pair in data:
        result = export_primer_pair(primer_pair)

    return result


def create_libamp_primer(pair, strand = "left") -> LibampPrimer:
    return LibampPrimer(
            pair[strand]["seq"],
            pair[strand]["gc_content"],
            pair[strand]["chr_start"],
            pair[strand]["chr_end"],
            pair[strand]["melting_temp"],
            strand,
            pair["score"],
            pair["product_size"],
            pair["version"],
            pair["targeton"],
            generate_primer_name(pair["pair"], strand),
        )


def generate_primer_name(pair_name, strand) -> str:
    append = "F" if strand == left else "R"

    return pair_name + append

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
                "value": ids["forward_direction"] if primer.strand == "left" else ids["negative_direction"],
            },
            "LibAmp Primer Type": {
                "value":  ids["libamp_forward"] if primer.strand == "left" else ids["libamp_reverse"],
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
        "folderId": ids["folder_id"],
        "schemaId": ids["libamp_schema_id"]
    }

def export_primer_pair(pair) -> None:
    benchling_ids = json.load(open("benchling_ids.json"))
    api_caller = Caller(benchling_connection.oligos_url)
    token = benchling_connection.token

    primer_left = create_libamp_primer(pair, "left")
    primer_right = create_libamp_primer(pair, "right")

    primer_left_json = primer_to_benchling_json(primer_left, benchling_ids)
    primer_right_json = primer_to_benchling_json(primer_right, benchling_ids)

    try:
        left_result = export_to_benchling(api_caller, token, primer_left_json)
        right_result = export_to_benchling(api_caller, token, primer_right_json)

        print("Result: ", left_result, right_result)

    except Exception as err:
        raise Exception(err)


def export_to_benchling(caller, token, json) -> str:
    res = caller.make_request('post', token, json).json()

    return res





