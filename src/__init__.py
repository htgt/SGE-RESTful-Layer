import sys
import os
from src.utils.exceptions import NoDotENVFile
from dotenv import load_dotenv

# Add src/.. to python path to find all modules inside.
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

try:
    load_dotenv(".env")
    BENCHLING_SECRET_KEY = os.getenv('BENCHLING_SECRET_KEY')
    BENCHLING_TENANT = os.getenv('BENCHLING_TENANT')
    GUNICORN_ENV = os.getenv('GUNICORN_ENV')
except:
    raise NoDotENVFile(f"No or invalid .env")
