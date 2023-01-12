import csv
import requests

class CSVReader():

    def parse_lines(self, url):
        content = self.load(url)
        lines = self.read(content)

        return lines

    def read(self, content):
        reader = csv.DictReader(content.splitlines(), delimiter=',')

        lines = []

        for row in reader:
            lines.append(row)

        return lines

    def load(self, url):
        s = requests.get(url).content.decode('utf-8')

        return s

