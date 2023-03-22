import unittest
from mock import MagicMock

from requests.models import Response

from src.benchling.push_libamp import export_primer_pair
from src.biology.libamp_primers import LibampPrimer


class ExportLibampPrimerPairToBenchling(unittest.TestCase):
    def setUp(self):
        self.left_primer = LibampPrimer(
            sequence='CTGTTCTGACAGTAGAAAGGCA',
            gc_content='45.45454545454545',
            chr_start='55',
            chr_end='77',
            melting_temp='58.004800503683725',
            strand='left',
            score='0.0',
            product_size=210,
            version='sfso_QEA2JL90',
            targeton='exon1',
            name='exon1_2_LibAmp_0F'
        )
        self.right_primer = LibampPrimer(
            sequence='AAGAATTTTCCCCAATGGTTGCT',
            gc_content='39.130434782608695',
            chr_start='242',
            chr_end='265',
            melting_temp='59.347613464584356',
            strand='right',
            score='0.0',
            product_size=210,
            version='sfso_QEA2JL90',
            targeton='exon1',
            name='exon1_2_LibAmp_0R'
        )

        self.fake_json = MagicMock()
        self.fake_json.return_value = {"id": "1"}

        self.good_response = MagicMock()
        self.good_response.status_code = 200
        self.good_response.ok = True
        self.good_response.json = self.fake_json

        self.bad_response = MagicMock()
        self.bad_response.status_code = 400
        self.bad_response.ok = False

    def test_left_primer_archive_if_right_fails(self):
        export_mock = MagicMock()
        export_mock.side_effect = [self.good_response, self.bad_response]

        archive_mock = MagicMock()

        test_url = 'http://www'
        test_token = 'abc'

        with self.assertRaises(Exception):
            export_primer_pair(
                self.left_primer,
                self.right_primer,
                export_mock,
                archive_mock,
                test_url,
                test_token
            )

        archive_mock.assert_called()

    def test_not_push_right_if_left_fails(self):
        archive_mock = MagicMock()

        export_mock = MagicMock()
        export_mock.return_value = self.bad_response

        test_url = 'http://www'
        test_token = 'abc'

        with self.assertRaises(Exception):
            export_primer_pair(
                self.left_primer,
                self.right_primer,
                export_mock,
                archive_mock,
                test_url,
                test_token
            )

        export_mock.called_only_once()
