from Bio.Seq import Seq
from dataclasses import dataclass
from src.utils.base_classes import BaseClass
from typing import List

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
        self.targeton = data['targeton']
        self.strand = data['strand']
        self.wge_link = data['wge_link']
        self.off_targets = data['off_targets']
        self.species = data['species']

    def as_benchling_req_body(self, event) -> dict:
        body = {
            'bases' : self.sequence,
            'fields': {
                'WGE ID' : {
                    'value' : self.wge_id,
                },
                'Guide Sequence' : {
                    'value' : self.sequence,
                },
                'Targeton' : {
                    'value' : self.targeton,
                },
                'Strand' : {
                    'value' : self.strand,
                },
                'WGE Hyperlink' : {
                    'value' : self.wge_link,
                },
                'Off Target Summary Data' : {
                    'value' : self.off_targets,
                },
                'Species' : {
                    'value' : self.species,
                },
            },
            'folderId' : event['folder_id'],
            'name' : event['name'],
            'schemaId' : event['schema_id'],
        }
        return body

@dataclass
class Oligo(BaseClass):
    sequence: str


@dataclass
class OligosPair(BaseClass):
    forward: Oligo
    reverse: Oligo
    
    def to_list_dicts(self) -> List[dict]:
        list_of_dicts = []
        for field in self.get_fields():
            return_dict = getattr(self, field)._asdict()
            list_of_dicts.append(return_dict)
        return list_of_dicts


class GuideRNAOligo(BaseClass):
    def __init__(self, seq) -> None:
        self.sequence = Seq(seq)
        self.first_base = transformations_dict["FIRST_BASE"]
        self.forward_prefix = transformations_dict["FORWARD_PREFIX"]
        self.reverse_prefix = transformations_dict["REVERSE_PREFIX"]
    
    def forward_sequence(self) -> Seq:
        return self.first_base + self.sequence[1:]
    
    def reverse_sequence(self) -> Seq:
        return self.forward_sequence().reverse_complement()

    def create_oligos(self) -> OligosPair:
        forward_oligo = Oligo(
            self.forward_prefix + self.forward_sequence(),
        )
        reverse_oligo = Oligo(
            self.reverse_prefix + self.reverse_sequence(),
        )

        return OligosPair(forward_oligo, reverse_oligo)
