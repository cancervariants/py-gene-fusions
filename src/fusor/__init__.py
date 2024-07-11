"""Provide computable object representation and validation for gene fusions"""

from importlib.metadata import PackageNotFoundError, version

from fusor.fusor import FUSOR

__all__ = ["__version__", "FUSOR"]


try:
    __version__ = version("cool_seq_tool")
except PackageNotFoundError:
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError
