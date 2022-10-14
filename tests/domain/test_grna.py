import unittest

from src.domain.guideRNA import GuideRNA

class TestGuideRNA(unittest.TestCase):
    def test_create_guide_RNA(self):
        input_data = {'id': '1168686327', 'sequence': 'GACTTCCAGCTACGGCGCG','gene_name': 'A1BG'}

        test_gRNA = GuideRNA(input_data)

        self.assertEqual(getattr(test_gRNA, "id"), "1168686327")
        self.assertEqual(getattr(test_gRNA, "sequence"), "GACTTCCAGCTACGGCGCG")
        self.assertEqual(getattr(test_gRNA, "gene_name"), "A1BG")

    #def test_create_set_of_guide_RNAs(self):

if __name__ == '__main__':
    unittest.main()
    