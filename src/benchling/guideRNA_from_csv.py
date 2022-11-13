from src.utils.file_system import CSVReader
from src.domain.guideRNA import create_set_of_gRNAs
from src.benchling.create_gRNA import export_grna_to_benchling

class GrnasImportFromCSV:
    def _get_lines_from_csv(self, url):
        lines = CSVReader().parse_lines(url)
        grnas_list = create_set_of_gRNAs(lines)

        return grnas_list

    def _post_to_benchling(self, items):
        for item in items:
            export_grna_to_benchling(item)

    def import_grnas(self, csv_url):
        grnas_list = self._get_lines_from_csv(csv_url)

        self._post_to_benchling(grnas_list)

        return grnas_list

