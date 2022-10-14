class GuideRNA:
    def __init__(self, data) -> None :
        setattr(self, "id", data["id"])
        setattr(self, "sequence", data["sequence"])
        setattr(self, "gene_name", data["gene_name"])
