from src.domain.taskImportGrnas import TaskImport
from src.rest_calls.send_calls import Caller
from src.benchling.guideRNA_from_csv import GrnasImportFromCSV

from src.benchling import benchling_connection

statuses = {
    "pending": "ttt",
    "failed_id": "ggg",
    "completed_id": "hhh",
}
TASKS_API_URL = 'https://tol-sangertest.benchling.com/api/v2/workflow-tasks/'
TASKS_OUTPUT_API_URL = 'https://tol-sangertest.benchling.com/api/v2/workflow-outputs'


class WorkflowTaskImport(TaskImport):
    def _get_task_update_url(self):
        return TASKS_API_URL + self.id

    def _get_status_id(self, id):
        return statuses[id]

    def execute(self):
        try:
            importer = GrnasImportFromCSV()
            result = importer.import_grnas(self.file_url)
        except Exception as err:
            raise Exception("Could not import guide RNAs")

        return result

    def update_status(self, status):
        url = self._get_task_update_url()

        print('URL:::::::::', url)

        api_caller = Caller(url)
        token = benchling_connection.token

        task_id = api_caller.make_request('patch', token, {"status": "wfts_DrqRAcOl"}).json()

        return task_id

    def add_task_output(self, payload):
        url = TASKS_OUTPUT_API_URL

        api_caller = Caller(url)
        token = benchling_connection.token

        output = self._prepare_task_output(self.id, payload)

        result = api_caller.make_request('post', token, output).json()

        return result

    def _prepare_task_output(self, task_id, grnas_list):
        json = {
            "fields": {
                "Oligos": {
                    "value": grnas_list
                },
            },
            "workflowTaskId": task_id
        }

        print("JSON ****", json)

        return json