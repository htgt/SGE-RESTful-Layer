from src.rest_calls.send_calls import Caller
from urllib.parse import urljoin

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
    def make_get(self, headers: dict, get_path: str) -> MockResponse: 
        url = urljoin(self.endpoint, get_path)
        return MockRequests.get(url, headers=headers)

    def make_post(self, headers: dict, json: dict) -> MockResponse: 
        return MockRequests.post(self.endpoint, json=json, headers=headers)
    
    def make_patch(self, headers: dict, json: dict) -> MockResponse: 
        return MockRequests.patch(self.endpoint, json=json, headers=headers)

def Mock_request_to_service(    
    json: dict,
    url: str,
    token: str,
    action: str = 'get'
) -> MockResponse:
    
    api_caller = MockCaller(url)
    response = api_caller.make_request(action, token, json)
    
    return response

def Mock_request_to_benchling(    
    service_url : str,
    action : str ,
    json = {}
) -> MockResponse:
    response = Mock_request_to_service(json, service_url, 'token', action=action)

    return response