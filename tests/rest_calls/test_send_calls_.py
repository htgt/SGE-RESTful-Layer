import unittest
import requests
from benchling import benchling_connection
from mock import patch

from rest_calls.send_calls import Caller


class TestCaller(unittest.TestCase):
    auth_mock = 'mocked_token'
    headers = {'Authorization': f"Bearer {auth_mock}"}

    def setUp(self):
        return

    @patch('requests.get')
    def test_make_get_success(self, request):
        # arrange
        request.return_value = requests.Response
        request.return_value.status_code = 404
        request.return_value.ok = True

        expected = requests.Response
        expected_header = {'Authorization': f"Bearer mocked_token"}

        test_endpoint = benchling_connection.api_url
        test_path = 'blobs/51cc7076-633d-42fc-a216-982fdc63a3ce'
        caller = Caller(test_endpoint)

        # act
        actual = caller.make_get(self.headers, test_path)

        # assert
        self.assertTrue(request.called)
        self.assertEqual(f"{request.call_args}",
                         f"call('{benchling_connection.blobs_url}51cc7076-633d-42fc-a216-982fdc63a3ce', headers={expected_header})")

    @patch('builtins.print')
    @patch('requests.post')
    def test_make_post_success(self, request, mock_print):
        # arrange
        request.return_value = requests.Response
        request.return_value.status_code = 200
        request.return_value.ok = True

        test_endpoint = 'endpoint'
        test_data = {'data_key': 'data_value'}
        expected_header = {'Authorization': f"Bearer mocked_token"}
        caller = Caller(test_endpoint)
        expected = requests.Response

        # act
        actual = caller.make_post(self.headers, test_data)

        # assert
        self.assertEqual(actual, expected)
        self.assertTrue(request.called)
        self.assertEqual(
            f"{request.call_args}",
            f"call('{test_endpoint}', json={test_data}, headers={expected_header})"
        )

    @patch('builtins.print')
    @patch('requests.post')
    def test_make_post_fail(self, request, mock_print):
        # arrange
        request.return_value = requests.Response
        request.return_value.status_code = 404
        request.return_value.reason = "Not Found"
        request.return_value.ok = False

        test_endpoint = 'end_point'
        test_data = {'data_key': 'data_value'}
        expected_header = {'Authorization': f"Bearer mocked_token"}
        caller = Caller(test_endpoint)
        expected = requests.Response

        # act
        actual = caller.make_post(self.headers, test_data)

        # assert
        self.assertEqual(actual, expected)
        self.assertTrue(mock_print.called)
        self.assertTrue(request.called)
        self.assertEqual(f"{request.call_args}",
                         f"call('{test_endpoint}', json={test_data}, headers={expected_header})")


if __name__ == '__main__':
    unittest.main()
