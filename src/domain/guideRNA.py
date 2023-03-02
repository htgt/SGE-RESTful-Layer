from Bio.Seq import Seq
from dataclasses import dataclass
from typing import List
from src.benchling.classes import BaseClass

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
class Oligo(BaseClass):
    sequence: str


@dataclass
class OligosPair(BaseClass):
    forward: Oligo
    reverse: Oligo

    def to_list_dicts(self) -> List[str]:
        list_of_dicts = []
        for oligo in self.get_fields():
            return_dict = vars(getattr(self, oligo))
            list_of_dicts.append(return_dict)
        return list_of_dicts


class GuideRNAOligo(BaseClass):
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
            self.forward_prefix + transformed_seq,
        )
        reverse_oligo = Oligo(
            self.reverse_prefix + transformed_seq.reverse_complement(),
        )

        return OligosPair(forward_oligo, reverse_oligo)
