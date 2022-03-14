"""Module containing methods and fixtures used throughout tests."""
from fusor import APP_ROOT
import pytest
import json

from fusor.fusor import FUSOR

EXAMPLES_DIR = APP_ROOT.resolve().parents[0] / "examples"


@pytest.fixture(scope="session")
def fusor():
    """Create test fixture for fusor object"""
    return FUSOR()


@pytest.fixture(scope="module")
def exhaustive_example():
    """Create test fixture for example of fusion (additional fields incl)"""
    return json.load(open(EXAMPLES_DIR / "exhaustive_example.json", "r"))


@pytest.fixture(scope="function")
def fusion_example():
    """Create test fixture for example of fusion (additional fields excl)"""
    return json.load(open(EXAMPLES_DIR / "minimal_example.json", "r"))
