import sys
import os
from src.utils.exceptions import MissingEnvVariables
from dotenv import load_dotenv
from pathlib import Path

# Add src/.. to python path to find all modules inside.
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

try:
    if Path(".env").is_file():
        print(".env file found, loading...")
        load_dotenv(".env")
    BENCHLING_SECRET_KEY = os.getenv('BENCHLING_SECRET_KEY')
    BENCHLING_TENANT = os.getenv('BENCHLING_TENANT')
    GUNICORN_ENV = os.getenv('GUNICORN_ENV')
    print(f"Benchling tenant: {BENCHLING_TENANT}")
except:
    raise MissingEnvVariables(f"No or invalid .env")

if 'unittest' in sys.modules or BENCHLING_TENANT == 'unittest':
    print("Unittest mode, changing tenant to 'unittest'")
    BENCHLING_SECRET_KEY = 'unittest'
    BENCHLING_TENANT = 'unittest'
    GUNICORN_ENV = 'unittest'
    