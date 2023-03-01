import unittest

from Bio.Seq import Seq
from src.domain.guideRNA import GuideRNA, GuideRNAOligo, create_set_of_gRNAs


class TestGuideRNA(unittest.TestCase):
    def test_create_guide_RNA(self):
        input_data = {'wge_id': '1168686327', 'seq': 'GACTTCCAGCTACGGCGCG', 'gene_symbol': 'A1BG'}

        test_gRNA = GuideRNA(input_data)

        self.assertEqual(getattr(test_gRNA, "id"), "1168686327")
        self.assertEqual(getattr(test_gRNA, "sequence"),
                         Seq("GACTTCCAGCTACGGCGCG"))
        self.assertEqual(getattr(test_gRNA, "gene_name"), "A1BG")

    def test_create_set_of_guide_RNAs(self):
        input_data = [{
            'wge_id': '1168686327',
            'seq': 'GACTTCCAGCTACGGCGCG',
            'gene_symbol'   : 'A1BG'
        }, {
            'wge_id'   : '1067960606',
            'seq': 'AATATGGTGGCCCTCCACC',
            'gene_symbol': 'A1CF'
        }]

        test_array = create_set_of_gRNAs(input_data)

        self.assertEqual(getattr(test_array[0], "id"), "1168686327")
        self.assertEqual(getattr(test_array[1], "id"), "1067960606")

    def test_forward_single_gRNA(self):
        input_data = {
            'wge_id'   : '1067960606',
            'seq': 'AATATGGTGGCCCTCCACC',
            'gene_symbol': 'A1CF'
        }

        test_gRNA = GuideRNA(input_data)

        self.assertEqual(test_gRNA.forward_sgRNA(), Seq("CACCAATATGGTGGCCCTCCACC"))

    def test_reverse_single_gRNA(self):
        input_data = {
            'wge_id': '1168686327',
            'seq': 'GACTTCCAGCTACGGCGCG',
            'gene_symbol'    : 'A1BG'
        }

        test_gRNA = GuideRNA(input_data)

        self.assertEqual(test_gRNA.reverse_sgRNA(), Seq("AAACCGCGCCGTAGCTGGAAGTC"))


class TestGuideRNAOligo(unittest.TestCase):
    def test_transformed_sequence_first_base(self):
        input_sequence = 'CACCAATATGGTGGCCCTCCATT'
        first_base = 'G'

        transformed = GuideRNAOligo(input_sequence).transform_first_and_last_bases()

        self.assertEqual(transformed[0], first_base)

    def test_transformed_sequence_last_base(self):
        input_sequence = 'CACCAATATGGTGGCCCTCCATT'
        last_base = 'C'

        transformed = GuideRNAOligo(input_sequence).transform_first_and_last_bases()

        self.assertEqual(transformed[-1], last_base)

    def test_create_oligos(self):
        input_sequence = 'AATATGGTGGCCCTCCATT'

        oligos = GuideRNAOligo(input_sequence).create_oligos()

        self.assertEqual(oligos.forward.sequence, Seq('CACCGATATGGTGGCCCTCCATC'))
        self.assertEqual(oligos.reverse.sequence, Seq('AAACGATGGAGGGCCACCATATC'))


if __name__ == '__main__':
    unittest.main()
