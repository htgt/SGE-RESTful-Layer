from flask import request
from flask_restful import Resource

import json

from src.benchling.workflow_task import WorkflowTaskImport

BENCHLING_WORKFLOW_TASK_SCHEMA_ID = "prstsch_bbPNDswA"
BENCHLING_UPDATED_STATUS_EVENT = "v2.workflowTask.updated.status"

class TaskEndpoint(Resource):
    def __transform_event_input_data(data):
        task_data = {}

        task_data["id"] = data["detail"]["workflowTask"]["id"]
        task_data["file_id"] = \
        data["detail"]["workflowTask"]["fields"]["gRNAs list CSV"]["value"]
        task_data["status_id"] = data["detail"]["workflowTask"]["status"]["id"]

        return task_data

    def get(self, id):

        return id, 201

    def post(self):
        data = request.json
        
        print("Task data received:")
        print(data)

        if data["detail-type"] == BENCHLING_UPDATED_STATUS_EVENT and data["detail"]["workflowTask"]["schema"]["id"] == BENCHLING_WORKFLOW_TASK_SCHEMA_ID :
            try:
                task_data = self.__transform_event_input_data(data)
                import_task = WorkflowTaskImport(task_data)

                created_grnas = import_task.execute()
                import_task.add_task_output(created_grnas)

                result = import_task.complete_task()

                return result, 200

            except Exception as err:
                return json.dumps(err), 500
        else:
            return "Incorrect input data", 404
