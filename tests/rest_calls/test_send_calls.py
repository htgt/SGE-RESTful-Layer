import unittest
import requests

from mock import patch

from rest_calls.send_calls import Caller


class TestCaller(unittest.TestCase):
    def setUp(self):
        return

    @staticmethod
    def auth_mock():
        return 'mocked_token'

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
        actual = caller.make_post(self.auth_mock, test_data)

        # assert
        self.assertEqual(actual, expected)
        self.assertTrue(request.called)
        self.assertEqual(f"{request.call_args}", f"call('{test_endpoint}', json={test_data}, headers={expected_header})")

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
        actual = caller.make_post(self.auth_mock, test_data)

        # assert
        self.assertEqual(actual, expected)
        self.assertTrue(mock_print.called)
        self.assertTrue(request.called)
        self.assertEqual(f"{request.call_args}", f"call('{test_endpoint}', json={test_data}, headers={expected_header})")

if __name__ == '__main__':
    unittest.main()
