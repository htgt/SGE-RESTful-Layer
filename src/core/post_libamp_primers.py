from src.benchling.push_libamp import export_primer_pair
from src.domain.libamp_primers import LibampPrimer


def post_libamp_primers(data):
    for primer_pair in data:
        primer_left = create_libamp_primer(primer_pair, "left")
        primer_right = create_libamp_primer(primer_pair, "right")

        result = export_primer_pair(primer_left, primer_right)

    return result

def generate_primer_name(pair_name, strand) -> str:
    append = "F" if strand == "left" else "R"

    return pair_name + append

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
