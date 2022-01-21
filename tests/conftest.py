"""Module containing methods and fixtures used throughout tests."""
from fusor import APP_ROOT
import pytest
import json

EXAMPLES_DIR = APP_ROOT.resolve().parents[0] / "examples"


@pytest.fixture(scope="module")
def exhaustive_example():
    """Create test fixture for example of fusion (additional fields incl)"""
    return json.load(open(EXAMPLES_DIR / "exhaustive_example.json", "r"))


@pytest.fixture(scope="module")
def fusion_example():
    """Create test fixture for example of fusion (additional fields excl)"""
    return json.load(open(EXAMPLES_DIR / "minimal_example.json", "r"))
