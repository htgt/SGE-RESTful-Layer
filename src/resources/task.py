from flask import request
from flask_restful import Resource

import json

from src.benchling.workflow_task import WorkflowTaskImport

BENCHLING_WORKFLOW_TASK_SCHEMA_ID = "prstsch_bbPNDswA"

class TaskEndpoint(Resource):
    def get(self, id):

        return id, 201


    def post(self):
       # data = request.json

        data = json.loads("""{
            "version": "0",
            "id": "70d32a26-20b2-81ae-fe94-5399384d269c",
            "detail-type": "v2.workflowTask.updated.status",
            "source": "aws.partner/benchling.com/tol-sangertest/sge-test",
            "account": "110091099050",
            "time": "2022-11-25T11:37:19Z",
            "region": "us-east-1",
            "resources": [],
            "detail": {
                "schema": {
                "id": "prstsch_bbPNDswA",
                "name": "SGE gRNA Auto Import"
            },
            "updates": [
                "status"
            ],
            "deprecated": false,
            "workflowTask": {
                "id": "wftask_qOtMxQl3",
                "fields": {
                    "gRNAs list CSV": {
                        "type": "blob_link",
                        "value": "747881dc-85b6-477f-97d4-c0acdb56bcde",
                        "isMulti": false,
                        "textValue": "gRNA_example.csv",
                        "displayValue": "gRNA_example.csv"
                    }
                },
                "schema": {
                    "id": "prstsch_bbPNDswA",
                    "name": "SGE gRNA Auto Import"
                },
                "status": {
                    "id": "wfts_EOjUQSei",
                    "statusType": "IN_PROGRESS",
                    "displayName": "In Progress"
                },
                "creator": {
                    "id": "ent_wsSRluG6",
                    "name": "Benchling Workflows",
                    "handle": ""
                },
                "assignee": null,
                "displayId": "sge-grna-auto2-T1",
                "clonedFrom": null,
                "scheduledOn": null,
                "archiveRecord": null,
                "executionType": "DIRECT",
                "creationOrigin": {
                    "originId": null,
                    "originType": "PROCESSES_SCHEMA_FIELD_MAPPER",
                    "application": "PROCESSES",
                    "originModalUuid": null
                },
                "executionOrigin": null,
                "workflowOutputs": [],
                "workflowTaskGroup": {
                    "id": "prs_97UJYcB0",
                    "name": "SGE gRNA Auto Import 2",
                    "displayId": "sge-grna-auto2"
                }
            },
            "excludedProperties": [],
            "id": "evt_aZriu5M5JM31",
            "createdAt": "2022-11-25T11:37:19.157644+00:00",
            "eventType": "v2.workflowTask.updated.status"
            }
        }""")

        def transform_data(json_data):
            task_data = {}
            task_data["id"] = data["detail"]["workflowTask"]["id"]
            task_data["file_id"] = data["detail"]["workflowTask"]["fields"]["gRNAs list CSV"]["value"]
            task_data["status_id"] = data["detail"]["workflowTask"]["status"]["id"]

            return task_data

        if data["detail-type"] == "v2.workflowTask.updated.status" and data["detail"]["workflowTask"]["schema"]["id"] == BENCHLING_WORKFLOW_TASK_SCHEMA_ID :

            try:
                task_data = transform_data(data)
                import_task = WorkflowTaskImport(task_data)

                created_grnas = import_task.execute()
                import_task.add_task_output(created_grnas)

                result = import_task.complete_task()

                return result, 200

            except Exception as err:
                return json.dumps(err), 500
        else:
            return "Incorrect input data", 404
