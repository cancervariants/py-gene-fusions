# FUSOR

FUSOR (**FUS**ion **O**bject **R**epresentation) provides modeling and validation tools for representing gene fusions in a flexible, computable structure.

### Installation

To install FUSOR:
```commandline
pip install fusor
```

For a development install, we recommend using Pipenv. See the
[pipenv docs](https://pipenv-fork.readthedocs.io/en/latest/#install-pipenv-today)
for direction on installing pipenv in your compute environment.

Once installed, from the project root dir, just run:

```commandline
pipenv lock
pipenv sync
```

### Data Downloads

#### SeqRepo
`FUSOR` relies on [seqrepo](https://github.com/biocommons/biocommons.seqrepo), which you must download yourself.

From the _root_ directory:
```
pip install seqrepo
sudo mkdir /usr/local/share/seqrepo
sudo chown $USER /usr/local/share/seqrepo
seqrepo pull -i 2021-01-29
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
