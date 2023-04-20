import unittest
from unittest.mock import patch
from src.benchling.connection.auth_utils import APIConnector
from tests.mock_connection import MockResponse


class TestAuthUtils(unittest.TestCase):
    def setUp(self):
        self.example = "example"
        self.example_token = "token"
        self.example_token_url = "http://test.com"
        self.example_client_id = "ID"
        self.example_key = "key"
        self.example_apiconnector_dict = {
            'token_url': self.example_token_url,
            'key': self.example_key,
            'client_id': self.example_client_id, 
            'auth_data': {'client_secret': self.example_key, 'client_id': self.example_client_id, 'grant_type': 'client_credentials'},
            'token': self.example_token
        }
    
    @patch('src.benchling.connection.auth_utils.APIConnector.get_access_token')
    def test_APIConnector(self, mocked_get_access_token):
        # Arrange
        expected_data = self.example_apiconnector_dict
        mocked_get_access_token.return_value = self.example_token
        # Act
        acted_data = APIConnector(self.example_token_url, self.example_client_id, self.example_key)
        # Assert
        self.assertDictEqual(vars(acted_data), expected_data)
        assert mocked_get_access_token.called
        
    @patch('requests.post')
    def test_get_access_token(self, mocked_post):
        # Arrange
        expected_data = self.example_apiconnector_dict
        mocked_post.return_value = MockResponse(200,'post', {"access_token":self.example_token})
        # Act
        acted_data = APIConnector(self.example_token_url, self.example_client_id, self.example_key)
        # Assert
        self.assertDictEqual(vars(acted_data), expected_data)
        assert mocked_post.called


if __name__ == '__main__':
    unittest.main()