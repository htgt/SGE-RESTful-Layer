from src.rest_calls.send_calls import Caller
from urllib.parse import urljoin
from typing import Tuple

class MockResponse:
    def __init__(self, status_code: int, action: str, json:dict = {}, ok = True) -> None:
        self.json_data = json
        self.status_code = status_code
        self.text = "test_text"
        self.ok = ok
        self.action = action

    def json(self):
        return self.json_data
        
class MockRequests:
    @staticmethod
    def get(url: str, *args, headers = {}, **kwargs): 
        if url == 'http://test.com' and headers:
            return MockResponse(200,'get')
        else:
            return MockResponse(404,'get', ok = False)

    @staticmethod
    def post(url: str, *args, json = {}, headers = {}, **kwargs):  
        if url == 'http://test.com' and json and headers:
            return MockResponse(200,'post' , json)
        else:
            return MockResponse(404,'post', json, ok = False)
    
    @staticmethod
    def patch(url: str, *args, json = {}, headers = {}, **kwargs): 
        if url == 'http://test.com' and json and headers:
            return MockResponse(200,'patch', json)
        else:
            return MockResponse(404,'patch', json, ok = False)
        
class MockCaller(Caller):
    def make_get(self, headers: dict, data: Tuple[dict, str]) -> MockResponse: 
        url = urljoin(self.endpoint, data)
        return MockRequests.get(url, headers=headers)

    def make_post(self, headers: dict, data: Tuple[dict, str]) -> MockResponse: 
        return MockRequests.post(self.endpoint, json=data, headers=headers)
    
    def make_patch(self, headers: dict, data: Tuple[dict, str]) -> MockResponse: 
        return MockRequests.patch(self.endpoint, json=data, headers=headers)

def Mock_request_to_service(    
    url: str,
    token: str,
    action: str,
    data: Tuple[dict, str]
) -> MockResponse:
    
    api_caller = MockCaller(url)
    response = api_caller.make_request(action, token, data)
    
    return response

def Mock_request_to_benchling(    
    service_url : str,
    action : str ,
    data: Tuple[dict, str]
) -> MockResponse:
    response = Mock_request_to_service(service_url, 'token', action, data)

    return response