import unittest

from domain.guideRNA import GuideRNA
import benchling.create_gRNA


class TestCreategRNA(unittest.TestCase):
    def setUp(self):
        return

    def test_prepare_sgrna_json_positive_success(self):
        # arrange
        test_grna = GuideRNA({
            'wge_id': 'grna_name',
            'seq': 'CAGTAGACACATGGTATTG',
            'gene_symbol': 'A1CF'
        })

        test_strand = '+'
        test_ids = {
            'positive_strand' : 'positive_strand_id',
            'folder_id': 'folder_id',
            'sgrna_schema_id': 'sgrna_schema_id'
        }
        expected = {
            "bases": 'CACCCAGTAGACACATGGTATTG',
            "fields": {
                "Strand": {
                    "value": 'positive_strand_id'
                }
            },
            "folderId": 'folder_id',
            "name": 'fwd_grna_name',
            "schemaId": 'sgrna_schema_id'
        }

        # act
        actual = benchling.create_gRNA.prepare_sgrna_json(test_grna, test_strand, test_ids)

        # assert
        self.assertEqual(actual, expected)

    def test_prepare_sgrna_json_negative_success(self):
        # arrange
        test_grna = GuideRNA({
            'wge_id': 'grna_name',
            'seq': 'CAGTAGACACATGGTATTG',
            'gene_symbol': 'A1CF'
        })

        test_strand = '-'
        test_ids = {
            'negative_strand': 'negative_strand_id',
            'folder_id': 'folder_id',
            'sgrna_schema_id': 'sgrna_schema_id'
        }
        expected = {
            "bases": 'AAACCAATACCATGTGTCTACTG',
            "fields": {
                "Strand": {
                    "value": 'negative_strand_id'
                }
            },
            "folderId": 'folder_id',
            "name": 'rev_grna_name',
            "schemaId": 'sgrna_schema_id'
        }

        # act
        actual = benchling.create_gRNA.prepare_sgrna_json(test_grna, test_strand, test_ids)

        # assert
        self.assertEqual(actual, expected)

    def test_prepare_grna_json_valid_success(self):
        # arrange
        test_grna = GuideRNA({
            'wge_id': '1234',
            'seq': 'CAGTAGACACATGGTATTG',
            'gene_symbol': 'A1CF'
        })
        test_fwd_sgrna_id = 'fwd_id'
        test_rev_sgrna_id = 'rev_id'
        test_ids = {
            'folder_id': 'folder_id',
            'grna_schema_id': 'grna_schema_id'
        }

        expected = {
            "bases": 'CAGTAGACACATGGTATTG',
            "fields": {
                "Gene Name": {
                    "value": 'A1CF',
                },
                "WGE ID": {
                    "value": 1234,
                },
                "Forward sgRNA": {
                    "value": 'fwd_id',
                },
                "Reverse sgRNA": {
                    "value": 'rev_id',
                },
            },
            "folderId": 'folder_id',
            "name": '1234',
            "schemaId": 'grna_schema_id'
        }

        # act
        actual = benchling.create_gRNA.prepare_grna_json(
            test_grna, 
            test_fwd_sgrna_id, 
            test_rev_sgrna_id, 
            test_ids
        )

        # assert
        self.assertEqual(actual, expected)
