# py-gene-fusions

Base Package for Fusion Curation

### Installation
For a development install, we recommend using Pipenv. See the 
[pipenv docs](https://pipenv-fork.readthedocs.io/en/latest/#install-pipenv-today) 
for direction on installing pipenv in your compute environment.
 
Once installed, from the project root dir, just run:

```commandline
pipenv lock
pipenv sync
```

To download py-gene-fusions, enter:
```commandline
python3 -m pip install -i https://test.pypi.org/simple/ py-gene-fusions
```

### Init coding style tests

Code style is managed by [flake8](https://github.com/PyCQA/flake8) and checked prior to commit.

We use [pre-commit](https://pre-commit.com/#usage) to run conformance tests.

This ensures:

* Check code style
* Check for added large files
* Detect AWS Credentials
* Detect Private Key

Before first commit run:

```commandline
pre-commit install
```


### Running unit tests

Running unit tests is as easy as pytest.

```commandline
pipenv run pytest
```