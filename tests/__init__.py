import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../tests')))

from src.benchling import BenchlingSchemaIds

benchling_ids = BenchlingSchemaIds('tests/fixtures/example_benchling_schema_ids.json').ids
