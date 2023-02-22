from Bio.Seq import Seq
from dataclasses import dataclass

transformations_dict = {
    "FORWARD_PREFIX": "CACC",
    "REVERSE_PREFIX": "AAAC",
    "FIRST_BASE": "G",
    "LAST_BASE": "C",
}

def create_set_of_gRNAs(data):
    set_of_gRNAs = []

    for line in data:
        set_of_gRNAs.append(GuideRNA(line))

    return set_of_gRNAs


class GuideRNA:
    def __init__(self, data) -> None :
        self.id = data["wge_id"]
        self.sequence = Seq(data["seq"])
        self.gene_name = data["gene_symbol"]

        self.forward_prefix = transformations_dict["FORWARD_PREFIX"]
        self.reverse_prefix = transformations_dict["REVERSE_PREFIX"]

    # forward and reverse will not be used as a part of guideRNA class
    def forward_sgRNA(self) -> Seq:
        return Seq(self.forward_prefix + self.sequence)

    def reverse_sgRNA(self) -> Seq:
        return Seq(self.reverse_prefix + self.sequence.reverse_complement())


    
@dataclass
class Oligo:
    sequence: str
    direction: str
    bases: list
    
@dataclass
class OligosPair:
    forward: Oligo
    reverse: Oligo

class GuideRNAOligo:
    def __init__(self, seq) -> None:
        self.sequence = Seq(seq)
        self.first_base = transformations_dict["FIRST_BASE"]
        self.last_base = transformations_dict["LAST_BASE"]
        self.forward_prefix = transformations_dict["FORWARD_PREFIX"]
        self.reverse_prefix = transformations_dict["REVERSE_PREFIX"]

    def transform_first_and_last_bases(self) -> str:
        return self.first_base + self.sequence[1:-1] + self.last_base

    def create_oligos(self) -> OligosPair:
        transformed_seq = self.transform_first_and_last_bases()
        forward_oligo = Oligo(
            transformed_seq, 
            self.forward_prefix, 
            [self.first_base, self.last_base]
        )
        reverse_oligo = Oligo(
            transformed_seq.reverse_complement(), 
            self.reverse_prefix, 
            [self.last_base, self.first_base]
        )

        return OligosPair(forward_oligo, reverse_oligo)
