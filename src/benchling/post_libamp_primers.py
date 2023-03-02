from src.domain.libamp_primers import LibampPrimer
import json
import sys
sys.path.append("..")

benchling_ids = json.load(open('benchling_ids.json'))

def post_libamp_primers(data):
    pairs = []

    for primer_pair in data:
        libamp_primer_left = create_libamp_primer(primer_pair, "left")
        libamp_primer_right = create_libamp_primer(primer_pair, "right")
        pairs.append(
            transform_primer_to_benchling_json(libamp_primer_left, benchling_ids)
        )
        pairs.append(
            transform_primer_to_benchling_json(libamp_primer_right, benchling_ids)
        )

    return pairs.__repr__()


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
        )

def transform_primer_to_benchling_json(primer, ids) -> dict:
    return {
        "bases": primer.sequence,
        "fields": {
            "GC Content (%)": {
                "value": primer.gc_content,
            },
            "Genome Location": {
                "value": primer.chr_start,
            },
            "LibAmp Primer Type": {
                "value": primer.strand == "left" if ids["libamp_forward"] else ids["libamp_reverse"],
            },
            "Primer Score": {
                "value": primer.score,
            },
            "Product Size (bp)":{
                "value": primer.product_size,
            },
            "Targeton": {
                "value": primer.targeton,
            },
            "Tm (°C)": {
                "value": primer.melting_temp,
            },
            "Version": {
                "value": primer.version,
            },
        },
        "folderId": ids["folder_id"],
        "schemaId": ids["libamp_schema_id"]
    }

