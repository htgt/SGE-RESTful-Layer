import requests

class Caller:
    def make_post(self, endpoint, data_dict):

        res = requests.post(endpoint, json=data_dict)

        if res.ok:
            return res
        else:
            print(f'Unsuccessful request. Status code: {res.status_code}. Reason: {res.reason}')
            return res

