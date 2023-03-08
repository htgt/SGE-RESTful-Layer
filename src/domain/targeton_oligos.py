class TargetonOligo:
    def __init__(self, data) -> None:
        self.ref_chromosome = data['chromosome']
        self.ref_strand = data['strand']
        self.ref_start = data['start']
        self.ref_end = data['end']
        self.r2_start = data['r2_start']
        self.r2_end = data['r2_end']
        self.action_vector = data['action_vector']
        self.ext_vector = data['ext_vector']
        self.sgrna_vector = data['sgrna_vector']
        self.version = data['version']
