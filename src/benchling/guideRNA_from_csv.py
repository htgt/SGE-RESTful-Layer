from src.utils.file_system import CSVReader
from src.biology.guideRNA import create_set_of_gRNAs
from src.benchling.create_gRNA import export_grna_to_benchling


class GrnasImportFromCSV:
    def _get_grnas_from_csv(self, url):
        lines = CSVReader().parse_lines(url)
        grnas_list = create_set_of_gRNAs(lines)

        return grnas_list

    def _post_to_benchling(self, items):
        sgrna_ids_list = []
        for item in items:
            new_ids = export_grna_to_benchling(item)
            sgrna_ids_list = [*sgrna_ids_list, *new_ids]

        return sgrna_ids_list

    def import_grnas(self, csv_url):
        grnas_list = self._get_grnas_from_csv(csv_url)

        created_sgrna_ids = self._post_to_benchling(grnas_list)

        return created_sgrna_ids

    def get_grnas(self, csv_url):
        grnas = self._get_grnas_from_csv(csv_url)
        result = []

        for item in grnas:
            result.append(item.spacer.__str__())

        return result
