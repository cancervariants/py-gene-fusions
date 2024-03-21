# FUSOR

[![image](https://img.shields.io/pypi/v/fusor.svg)](https://pypi.python.org/pypi/fusor)
[![image](https://img.shields.io/pypi/l/fusor.svg)](https://pypi.python.org/pypi/fusor)
[![image](https://img.shields.io/pypi/pyversions/fusor.svg)](https://pypi.python.org/pypi/fusor)
[![Actions status](https://github.com/cancervariants/fusor/actions/workflows/checks.yaml/badge.svg)](https://github.com/cancervariants/fusor/actions/checks.yaml)

<!-- description -->
FUSOR (**FUS**ion **O**bject **R**epresentation) provides modeling and validation tools for representing gene fusions in a flexible, computable structure.
<!-- /description -->

## Installation

### Package

To install FUSOR:
```shell
python3 -m pip install fusor
```

### SeqRepo and UTA

FUSOR relies on data from [SeqRepo](https://github.com/biocommons/seqrepo) and the [Universal Transcript Archive (UTA)](https://github.com/biocommons/uta), supplied via [Cool-Seq-Tool](https://github.com/genomicmedlab/cool-seq-tool). See the [Cool-Seq-Tool installation instructions](https://coolseqtool.readthedocs.io/en/latest/install.html) for specifics on setting up SeqRepo and UTA.

To connect to the UTA database, you can use the default url (`postgresql://uta_admin@localhost:5433/uta/uta_20210129`). If you use the default url, you must either set the password using environment variable `UTA_PASSWORD` or setting the parameter `db_pwd` in the `FUSOR` class.

If you do not wish to use the default, you must set the environment variable `UTA_DB_URL` which has the format of `driver://user:pass@host:port/database/schema`.

### Gene Normalizer

`FUSOR` also relies on data from the [Gene Normalizer](https://github.com/cancervariants/gene-normalization), which you must download yourself. See the [README](https://github.com/cancervariants/gene-normalization#readme) for deploying the gene database.

## Development

### Setup

Clone the repo:

```shell
git clone https://github.com/cancervariants/fusor
cd fusor
```

Create a virtual environment, and activate it:

```shell
python3 -m virtualenv venv
source venv/bin/activate
```

Install test and dev dependencies:

```shell
python3 -m pip install -e '.[dev,tests]'
```

### Style

Code style is managed by [Ruff](https://github.com/astral-sh/ruff) and checked prior to commit.

```shell
python3 -m ruff format . && python3 -m ruff check --fix .
```

We use [pre-commit](https://pre-commit.com/#usage) to run conformance tests.

This performs checks for:

* Code style
* File endings
* Added large files
* AWS credentials
* Private keys

Before your first commit, run:

```shell
pre-commit install
```

### Unit tests

Unit testing is provided via `pytest`:

```shell
pytest
```
