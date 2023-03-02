from src.domain.libamp_primers import LibampPrimer


def post_libamp_primers(data):
    pairs = []

    for primer_pair in data:
        libamp_primer_left = create_libamp_primer(primer_pair, "left")
        libamp_primer_right = create_libamp_primer(primer_pair, "right")
        pairs.append(libamp_primer_left)
        pairs.append(libamp_primer_right)

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

