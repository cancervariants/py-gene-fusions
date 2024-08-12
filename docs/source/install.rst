.. _install:

Installation
============

Prerequisites
-------------

* Python 3.10+
* A UNIX-like environment (e.g. MacOS, WSL, Ubuntu)
* A recent version of PostgreSQL (ideally at least 11+)
* A modern Java runtime (if using DynamoDB for the Gene Normalizer database)

Library installation
--------------------

Install ``FUSOR`` from `PyPI <https://pypi.org/project/fusor/>`_:

.. code-block:: shell

    pip install fusor

Data setup
----------

Universal Transcript Archive (UTA)
++++++++++++++++++++++++++++++++++

The `UTA <https://github.com/biocommons/uta>`_ is a dataset of genome-transcript aligned data supplied as a PostgreSQL database. Access in FUSOR is supplied by way of ``Cool-Seq-Tool``; see the `Cool-Seq-Tool UTA docs <https://coolseqtool.readthedocs.io/stable/install.html#set-up-uta>`_ for some opinionated setup instructions.

SeqRepo
+++++++

`SeqRepo <https://github.com/biocommons/biocommons.seqrepo>`_ is a controlled dataset of biological sequences. As with UTA, access in FUSOR is given via `Cool-Seq-Tool`, which provides `documentation <https://coolseqtool.readthedocs.io/stable/install.html#set-up-seqrepo>`_ on getting it set up.

Gene Normalizer
+++++++++++++++

Finally, ``FUSOR`` uses the `Gene Normalizer <https://github.com/cancervariants/gene-normalization>`_ to ground gene terms. See the `Gene Normalizer documentation <https://gene-normalizer.readthedocs.io/stable/install.html>`_ for setup instructions.
