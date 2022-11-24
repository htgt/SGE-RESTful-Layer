from src.domain.taskImportGrnas import TaskImport
from src.rest_calls.send_calls import Caller

from src.benchling import benchling_connection

statuses = {
    "pending": "ttt",
    "failed_id": "ggg",
    "completed_id": "hhh",
}
tasks_api_url = 'https://tol-sangertest.benchling.com/api/v2/workflow-tasks/'


class WorkflowTaskImport(TaskImport):
    def _get_task_update_url(self):
        return tasks_api_url + self.id

    def _get_status_id(self, id):
        return statuses[id]

    def update_status(self, status):
        url = self._get_task_update_url()

        print('URL:::::::::', url)

        api_caller = Caller(url)
        token = benchling_connection.token

        task_id = api_caller.make_request('patch', token, {"status": "wfts_DrqRAcOl"}).json()

        return task_id