import unittest

from Bio.Seq import Seq
from src.domain.guideRNA import GuideRNA, create_set_of_gRNAs

class TestGuideRNA(unittest.TestCase):
    def test_create_guide_RNA(self):
        input_data = {'id': '1168686327', 'sequence': 'GACTTCCAGCTACGGCGCG','gene_name': 'A1BG'}

        test_gRNA = GuideRNA(input_data)

        self.assertEqual(getattr(test_gRNA, "id"), "1168686327")

    def test_create_set_of_guide_RNAs(self):
        input_data = [{
                'id': '1168686327',
                'sequence': 'GACTTCCAGCTACGGCGCG',
                'gene_name'   : 'A1BG'
            }, {
                'id'   : '1067960606',
                'sequence': 'AATATGGTGGCCCTCCACC',
                'gene_name': 'A1CF'
            }]

        test_array = create_set_of_gRNAs(input_data)

        self.assertEqual(getattr(test_array[0], "id"), "1168686327")
        self.assertEqual(getattr(test_array[1], "id"), "1067960606")

    def test_forward_single_gRNA(self):
        input_data = {
                'id'   : '1067960606',
                'sequence': 'AATATGGTGGCCCTCCACC',
                'gene_name': 'A1CF'
            }

        test_gRNA = GuideRNA(input_data)

        self.assertEqual(test_gRNA.forward_sgRNA(), Seq("CACCAATATGGTGGCCCTCCACC"))

    def test_reverse_single_gRNA(self):
        input_data = {
            'id': '1168686327',
            'sequence': 'GACTTCCAGCTACGGCGCG',
            'gene_name'    : 'A1BG'
        }

        test_gRNA = GuideRNA(input_data)

        self.assertEqual(test_gRNA.reverse_sgRNA(), Seq("AAACCGCGCCGTAGCTGGAAGTC"))

if __name__ == '__main__':
    unittest.main()
    