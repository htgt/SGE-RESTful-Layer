class MockResponse:
    def __init__(self, json_data, status_code, ok = True):
        self.json_data = json_data
        self.status_code = status_code
        self.text = "test_text"
        self.ok = ok

    def json(self):
        return self.json_data
        
class MockRequests:
    @staticmethod
    def get(json: dict = {}, *args, **kwargs): 
        if args[0] == 'http://test.com':
            return MockResponse(json, 200)
        else:
            return MockResponse(json, 404, ok = False)

    @staticmethod
    def post(json: dict = {}, *args, **kwargs): 
        if args[0] == 'http://test.com':
            return MockResponse(json, 200)
        else:
            return MockResponse(json, 404, ok = False)
    
    @staticmethod
    def patch(json: dict = {}, *args, **kwargs): 
        if args[0] == 'http://test.com':
            return MockResponse(json, 200)
        else:
            return MockResponse(json, 404, ok = False)
        
