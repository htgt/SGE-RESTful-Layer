from dataclasses import dataclass


@dataclass
class LibampPrimer:
    sequence: str
    gc_content: float
    chr_start: int
    chr_end: int
    melting_temp: float
    strand: str
    score: float
    product_size: float
    version: str
    targeton: str
    name: str




