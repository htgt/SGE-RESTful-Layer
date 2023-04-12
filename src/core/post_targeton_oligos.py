from __future__ import annotations
from typing import TYPE_CHECKING

from src.benchling.utils.export_to_benchling import export_to_benchling
import json

if TYPE_CHECKING:
    from src.benchling.connection.connection_class import BenchlingConnection

def send_targeton_oligo_post_request(body: dict, benchling_connection: BenchlingConnection) -> dict:
    try:
        response = export_to_benchling(
            body,
            benchling_connection.custom_entity_url,
            benchling_connection,
            'post'
        )
        return response, 201
    except Exception as err:
        print(json.dumps(err))
        return json.dumps(err), 500
