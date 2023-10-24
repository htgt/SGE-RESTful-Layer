from src.biology.guideRNA import GuideRNA
from src.benchling.utils.request_to_benchling import request_to_benchling_json_response
from src.benchling import benchling_schema_ids, benchling_urls


def patch_guide_rna(guide: GuideRNA, event_data: dict) -> str:
    url = benchling_urls.guide_rna_url
    benchling_body = as_benchling_req_body(guide, event_data)
    patch_url = url + '/' + event_data["entity_id"]
    
    response = request_to_benchling_json_response(
        patch_url,
        'patch',
        benchling_body
    )

    return response

def as_benchling_req_body(guide: GuideRNA, event: dict) -> dict:
    species_benchling_id = benchling_schema_ids.ids["dropdowns"]["species"][guide.species]

    body = {
        'bases': guide.spacer,
        'fields': {
            'WGE ID': {
                'value': guide.wge_id,
            },
            'Targeton': {
                'value': event['targeton_id'],
            },
            #'Strand': {
            #   'value' : self.strand,
            #},
            'WGE Hyperlink': {
                'value': guide.wge_link,
            },
            'Off Target Summary Data': {
                'value': guide.off_targets,
            },
            'Species': {
                'value': species_benchling_id,
            },
            'PAM Sequence': {
                'value': guide.pam
            },
        },
        'name': event['name'],
        'schemaId' : event['schema_id'],
    }

    return body
