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
    print(f"Benchling tenant: {BENCHLING_TENANT}")
except:
    raise MissingEnvVariables(f"No or invalid .env")