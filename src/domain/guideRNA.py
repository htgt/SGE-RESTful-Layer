def create_set_of_gRNAs(data):
    set_of_gRNAs = []

    for line in data:
        set_of_gRNAs.append(GuideRNA(line))

    return set_of_gRNAs

class GuideRNA:
    FORWARD_PREFIX = "CACC"
    REVERSE_PREFIX = "AAAC"

    def __init__(self, data) -> None :
        setattr(self, "id", data["id"])
        setattr(self, "sequence", data["sequence"])
        setattr(self, "gene_name", data["gene_name"])

    def forward_sgRNA(self) -> str:
        return self.FORWARD_PREFIX + getattr(self, "sequence")
