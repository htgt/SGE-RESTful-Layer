from Bio.Seq import Seq
from dataclasses import dataclass
from src.utils.base_classes import BaseClass
from src.domain.species import get_species_name_by_id


transformations_dict = {
    "FORWARD_PREFIX": "CACC",
    "REVERSE_PREFIX": "AAAC",
    "FIRST_BASE": "G",
}


def create_set_of_gRNAs(data):
    set_of_gRNAs = []

    for line in data:
        set_of_gRNAs.append(GuideRNA(line))

    return set_of_gRNAs


class GuideRNA:
    def __init__(self, data) -> None :
        self.wge_id = data['wge_id']
        self.sequence = data['seq']
        #self.targeton = data['targeton']
        #self.strand = data['strand']
        self.wge_link = data['wge_link']
        self.off_targets = data['off_targets']
        self.species = get_species_name_by_id(data['species'])


@dataclass
class Oligo(BaseClass):
    sequence: str


class GuideRNAOligos(BaseClass):
    def __init__(self, seq) -> None:
        self.sequence = Seq(seq)
        self.first_base = transformations_dict["FIRST_BASE"]
        self.forward_prefix = transformations_dict["FORWARD_PREFIX"]
        self.reverse_prefix = transformations_dict["REVERSE_PREFIX"]

        self.forward = Oligo(
            self.forward_prefix + self.forward_sequence(),
        )
        self.reverse = Oligo(
            self.reverse_prefix + self.reverse_sequence(),
        )
    
    def forward_sequence(self) -> Seq:
        return self.first_base + self.sequence[1:]
    
    def reverse_sequence(self) -> Seq:
        return self.forward_sequence().reverse_complement()

