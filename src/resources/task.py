from flask import request
from flask_restful import Resource

import json

from src.domain.taskImportGrnas import TaskImport

class TaskEndpoint(Resource):
    def get(self, id):
        return id, 201


    def post(self):
       # data = request.json

        data = json.loads("""{
            "id": "cda06ef4-eac1-9e84-8777-c8f9cad34605",
            "detail-type": "v2.workflowTask.updated.status",
            "detail": {
                "workflowTask": {
                    "id": "wftask_dMvTXfUn",
                    "fields": {
                        "gRNAs list CSV": {
                            "type": "blob_link",
                            "value": "773b18f1-1658-4865-aee3-f512a8d0e869"
                        }
                    },
                    "schema": {
                        "id": "prstsch_xlEsqy9T",
                        "name": "SGE Import gRNAs"
                    },
                    "status": {
                        "id": "wfts_DrqRAcOl",
                        "statusType": "PENDING",
                        "displayName": "Pending"
                    },
                    "workflowOutputs": [],
                    "workflowTaskGroup": {
                        "id": "prs_W2ENmFcb",
                        "name": "SGE Import gRNAs 9",
                        "displayId": "sge-grna9"
                    }
                },
                "id": "evt_JmYcUA0O29oq",
                "createdAt": "2022-11-07T13:19:46.375574+00:00",
                "eventType": "v2.workflowTask.updated.status"
            }
        }""")

        def transform_data(json_data):
            task_data = {}
            task_data["id"] = data["detail"]["workflowTask"]["id"]
            task_data["file_id"] = data["detail"]["workflowTask"]["fields"]["gRNAs list CSV"]["value"]
            task_data["status_id"] = data["detail"]["workflowTask"]["status"]["id"]

            return task_data

        if data["detail-type"] == "v2.workflowTask.updated.status" and data["detail"]["workflowTask"]["schema"]["id"] == "prstsch_xlEsqy9T" :
            task_data = transform_data(data)

            import_task = TaskImport(task_data)
            execute_result = import_task.execute()

            return execute_result
        else:
            return "Incorrect input data", 404