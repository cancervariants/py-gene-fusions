"""Fusion package"""
import logging
from os import environ
from pathlib import Path

from .version import __version__

logging.basicConfig(
    filename="fusor.log",
    format="[%(asctime)s] - %(name)s - %(levelname)s : %(message)s",
)
logger = logging.getLogger("fusor")
logger.setLevel(logging.DEBUG)


SEQREPO_DATA_PATH = environ.get("SEQREPO_DATA_PATH", "/usr/local/share/seqrepo/latest")

if "UTA_DB_URL" in environ:
    UTA_DB_URL = environ["UTA_DB_URL"]
else:
    UTA_DB_URL = "postgresql://uta_admin@localhost:5433/uta/uta_20210129"

logging.getLogger("boto3").setLevel(logging.INFO)
logging.getLogger("botocore").setLevel(logging.INFO)
logging.getLogger("nose").setLevel(logging.INFO)
logging.getLogger("python_jsonschema_objects.classbuilder").setLevel(logging.INFO)
logging.getLogger("urllib3").setLevel(logging.INFO)

from fusor.fusor import FUSOR  # noqa: E402

APP_ROOT = Path(__file__).resolve().parents[0]
