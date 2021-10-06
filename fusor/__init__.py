"""Fusion package"""
from pathlib import Path
from os import environ

if "SEQREPO_DATA_PATH" in environ:
    SEQREPO_DATA_PATH = environ["SEQREPO_DATA_PATH"]
else:
    SEQREPO_DATA_PATH = "/usr/local/share/seqrepo/latest"

from fusor.fusor import FUSOR  # noqa: E402, F401

APP_ROOT = Path(__file__).resolve().parents[0]
