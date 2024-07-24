"""Module containing methods and fixtures used throughout tests."""

import logging

import pytest
from cool_seq_tool.app import CoolSeqTool

from fusor.fusor import FUSOR


def pytest_addoption(parser):
    """Add custom commands to pytest invocation.
    See https://docs.pytest.org/en/7.1.x/reference/reference.html#parser
    """
    parser.addoption(
        "--verbose-logs",
        action="store_true",
        default=False,
        help="show noisy module logs",
    )


def pytest_configure(config):
    """Configure pytest setup."""
    if not config.getoption("--verbose-logs"):
        logging.getLogger("botocore").setLevel(logging.INFO)
        logging.getLogger("boto3").setLevel(logging.INFO)
        logging.getLogger("urllib3").setLevel(logging.INFO)
        logging.getLogger("nose").setLevel(logging.INFO)


@pytest.fixture(scope="session")
def fusor_instance():
    """Create test fixture for fusor object

    Suppresses checks for CoolSeqTool external resources. Otherwise, on CST startup,
    it will try to check that its MANE summary file is up-to-date, which is an FTP call
    to the NCBI servers and can hang sometimes.

    If those files aren't available, create a CST instance in another session -- by
    default, it should save files to a centralized location that this test instance can
    access.
    """
    cst = CoolSeqTool(force_local_files=True)
    return FUSOR(cool_seq_tool=cst)


@pytest.fixture(scope="session")
def braf_gene_descriptor():
    """Create gene descriptor params for BRAF."""
    return {
        "type": "Gene",
        "id": "normalize.gene.hgnc:1097",
        "label": "BRAF",
        "mappings": [
            {
                "coding": {"code": "673", "system": "ncbigene"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "ENSG00000157764", "system": "ensembl"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "CCDS5863", "system": "ccds"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "1943", "system": "iuphar"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "119066", "system": "orphanet"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "BRAF", "system": "cosmic"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "2284096", "system": "pubmed"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "uc003vwc.5", "system": "ucsc"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "164757", "system": "omim"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "NM_004333", "system": "refseq"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "CCDS87555", "system": "ccds"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "P15056", "system": "uniprot"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "M95712", "system": "ena.embl"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "OTTHUMG00000157457", "system": "vega"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "1565476", "system": "pubmed"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "CCDS94219", "system": "ccds"},
                "relation": "relatedMatch",
            },
            {
                "coding": {"code": "CCDS94218", "system": "ccds"},
                "relation": "relatedMatch",
            },
        ],
        "alternativeLabels": ["BRAF1", "BRAF-1", "RAFB1", "NS7", "B-RAF1", "B-raf"],
        "extensions": [
            {
                "name": "approved_name",
                "value": "B-Raf proto-oncogene, serine/threonine kinase",
            },
            {
                "name": "ensembl_locations",
                "value": [
                    {
                        "type": "SequenceLocation",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        },
                        "start": 140719326,
                        "end": 140924929,
                    }
                ],
            },
            {
                "name": "ncbi_locations",
                "value": [
                    {
                        "type": "SequenceLocation",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        },
                        "start": 140713327,
                        "end": 140924929,
                    }
                ],
            },
            {"name": "ncbi_gene_type", "value": "protein-coding"},
            {
                "name": "hgnc_locus_type",
                "value": "gene with protein product",
            },
            {"name": "ensembl_biotype", "value": "protein_coding"},
            {"name": "strand", "value": "-"},
            {"name": "symbol_status", "value": "approved"},
        ],
    }


@pytest.fixture(scope="session")
def alk_gene_descriptor():
    """Create test fixture for ALK gene descriptor params"""
    return {
        "id": "normalize.gene.hgnc:427",
        "type": "Gene",
        "label": "ALK",
        "description": None,
        "alternativeLabels": ["NBLST3", "CD246", "ALK1"],
        "extensions": [
            {"name": "symbol_status", "value": "approved", "description": None},
            {
                "name": "approved_name",
                "value": "ALK receptor tyrosine kinase",
                "description": None,
            },
            {"name": "strand", "value": "-", "description": None},
            {
                "name": "ensembl_locations",
                "value": [
                    {
                        "id": "ga4gh:SL.V-yTsF-F4eHxeDHeU5KZIF3ZOzE2vUnG",
                        "type": "SequenceLocation",
                        "digest": "V-yTsF-F4eHxeDHeU5KZIF3ZOzE2vUnG",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
                        },
                        "start": 29192773,
                        "end": 29921586,
                    }
                ],
                "description": None,
            },
            {
                "name": "ncbi_locations",
                "value": [
                    {
                        "id": "ga4gh:SL.V-yTsF-F4eHxeDHeU5KZIF3ZOzE2vUnG",
                        "type": "SequenceLocation",
                        "digest": "V-yTsF-F4eHxeDHeU5KZIF3ZOzE2vUnG",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
                        },
                        "start": 29192773,
                        "end": 29921586,
                    }
                ],
                "description": None,
            },
            {
                "name": "hgnc_locus_type",
                "value": "gene with protein product",
                "description": None,
            },
            {"name": "ncbi_gene_type", "value": "protein-coding", "description": None},
            {"name": "ensembl_biotype", "value": "protein_coding", "description": None},
        ],
        "mappings": [
            {
                "coding": {
                    "label": None,
                    "system": "ensembl",
                    "version": None,
                    "code": "ENSG00000171094",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ncbigene",
                    "version": None,
                    "code": "238",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "orphanet",
                    "version": None,
                    "code": "160020",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "hcdmdb",
                    "version": None,
                    "code": "CD246",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ucsc",
                    "version": None,
                    "code": "uc002rmy.4",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "refseq",
                    "version": None,
                    "code": "NM_004304",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ccds",
                    "version": None,
                    "code": "CCDS33172",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "omim",
                    "version": None,
                    "code": "105590",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ena.embl",
                    "version": None,
                    "code": "D45915",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "vega",
                    "version": None,
                    "code": "OTTHUMG00000152034",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "uniprot",
                    "version": None,
                    "code": "Q9UM73",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "iuphar",
                    "version": None,
                    "code": "1839",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "pubmed",
                    "version": None,
                    "code": "8122112",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "cosmic",
                    "version": None,
                    "code": "ALK",
                },
                "relation": "relatedMatch",
            },
            {
                "coding": {
                    "label": None,
                    "system": "ccds",
                    "version": None,
                    "code": "CCDS86828",
                },
                "relation": "relatedMatch",
            },
        ],
    }


@pytest.fixture(scope="module")
def exhaustive_example(alk_gene_descriptor, braf_gene_descriptor):
    """Create test fixture for a fake fusion exemplifying most major field types, in
    'expanded' form (ie properties augmented by VICC descriptors)
    """
    return {
        "type": "CategoricalFusion",
        "critical_functional_domains": [
            {
                "type": "FunctionalDomain",
                "_id": "interpro:IPR020635",
                "label": "Tyrosine-protein kinase, catalytic domain",
                "status": "lost",
                "associatedGene": alk_gene_descriptor,
                "sequence_location": {
                    "id": "fusor.location_descriptor:NP_004295.2",
                    "type": "LocationDescriptor",
                    "label": None,
                    "description": None,
                    "xrefs": None,
                    "alternate_labels": None,
                    "extensions": None,
                    "location_id": "ga4gh:VSL.hQKhk6ZOOYZAmShXrzhfb6H3j65ovsKu",
                    "location": {
                        "id": None,
                        "type": "SequenceLocation",
                        "sequence_id": "refseq:NP_004295.2",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 1116},
                            "end": {"type": "Number", "value": 1383},
                        },
                    },
                },
            }
        ],
        "structure": [
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_152263.3",
                "exon_start": 1,
                "exon_start_offset": 0,
                "exon_end": 8,
                "exon_end_offset": 0,
                "gene_descriptor": {
                    "id": "normalize.gene:TPM3",
                    "type": "Gene",
                    "label": "TPM3",
                    "description": None,
                    "xrefs": ["ensembl:ENSG00000143549", "ncbigene:7170"],
                    "alternate_labels": [
                        "TM-5",
                        "TM5",
                        "NEM1~withdrawn",
                        "OK/SW-cl.5",
                        "TM30nm",
                        "TPMsk3",
                        "HEL-S-82p",
                        "TRK",
                        "CAPM1",
                        "TPM3nu",
                        "FLJ35371",
                        "TM30",
                        "TM3",
                        "CFTD",
                        "NEM1",
                        "hscp30",
                        "HEL-189",
                    ],
                    "extensions": [
                        {
                            "type": "Extension",
                            "name": "symbol_status",
                            "value": "approved",
                        },
                        {
                            "type": "Extension",
                            "name": "approved_name",
                            "value": "tropomyosin 3",
                        },
                        {
                            "type": "Extension",
                            "name": "hgnc_locations",
                            "value": [
                                {
                                    "_id": "ga4gh:VCL.rmJvYV5JccRSEoMVxe5BmuHs9S2VZ4uR",
                                    "type": "ChromosomeLocation",
                                    "species_id": "taxonomy:9606",
                                    "chr": "1",
                                    "interval": {
                                        "end": "q21.3",
                                        "start": "q21.3",
                                        "type": "CytobandInterval",
                                    },
                                }
                            ],
                        },
                        {
                            "type": "Extension",
                            "name": "ensembl_locations",
                            "value": [
                                {
                                    "_id": "ga4gh:VSL._ASa2-iBSDZSpC3JlpwJxzv4OY5M-5Ct",
                                    "type": "SequenceLocation",
                                    "sequence_id": "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                                    "interval": {
                                        "start": {"type": "Number", "value": 154155307},
                                        "end": {"type": "Number", "value": 154194648},
                                        "type": "SequenceInterval",
                                    },
                                }
                            ],
                        },
                        {
                            "type": "Extension",
                            "name": "ncbi_locations",
                            "value": [
                                {
                                    "_id": "ga4gh:VCL.rmJvYV5JccRSEoMVxe5BmuHs9S2VZ4uR",
                                    "type": "ChromosomeLocation",
                                    "species_id": "taxonomy:9606",
                                    "chr": "1",
                                    "interval": {
                                        "end": "q21.3",
                                        "start": "q21.3",
                                        "type": "CytobandInterval",
                                    },
                                },
                                {
                                    "_id": "ga4gh:VSL.sGJqQhhTg3BYlndAP7nFzN7KoKID1yP_",
                                    "type": "SequenceLocation",
                                    "sequence_id": "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                                    "interval": {
                                        "start": {"type": "Number", "value": 154155307},
                                        "end": {"type": "Number", "value": 154192100},
                                        "type": "SequenceInterval",
                                    },
                                },
                            ],
                        },
                        {
                            "type": "Extension",
                            "name": "associated_with",
                            "value": [
                                "vega:OTTHUMG00000035853",
                                "ccds:CCDS41400",
                                "ccds:CCDS1060",
                                "ccds:CCDS41402",
                                "ccds:CCDS41401",
                                "pubmed:25369766",
                                "cosmic:TPM3",
                                "refseq:NM_152263",
                                "orphanet:120227",
                                "uniprot:P06753",
                                "ccds:CCDS72922",
                                "ccds:CCDS60274",
                                "ucsc:uc001fec.3",
                                "omim:191030",
                                "ccds:CCDS41403",
                                "ena.embl:BC008425",
                                "ccds:CCDS60275",
                                "pubmed:1829807",
                            ],
                        },
                        {
                            "type": "Extension",
                            "name": "previous_symbols",
                            "value": ["NEM1", "NEM1~withdrawn", "FLJ35371"],
                        },
                        {
                            "type": "Extension",
                            "name": "hgnc_locus_type",
                            "value": "gene with protein product",
                        },
                        {
                            "type": "Extension",
                            "name": "ncbi_gene_type",
                            "value": "protein-coding",
                        },
                        {
                            "type": "Extension",
                            "name": "ensembl_biotype",
                            "value": "protein_coding",
                        },
                        {"type": "Extension", "name": "strand", "value": "-"},
                    ],
                    "gene_id": "hgnc:12012",
                    "gene": None,
                },
                "elementGenomicStart": {
                    "id": "fusor.location_descriptor:NC_000001.11",
                    "type": "LocationDescriptor",
                    "label": "NC_000001.11",
                    "description": None,
                    "xrefs": None,
                    "alternate_labels": None,
                    "extensions": None,
                    "location_id": "ga4gh:VSL.n7i6VMRAuSgAjwVopxhWAJdlPJMfk7KR",
                    "location": {
                        "id": None,
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 154192135},
                            "end": {"type": "Number", "value": 154192136},
                        },
                    },
                },
                "elementGenomicEnd": {
                    "id": "fusor.location_descriptor:NC_000001.11",
                    "type": "LocationDescriptor",
                    "label": "NC_000001.11",
                    "description": None,
                    "xrefs": None,
                    "alternate_labels": None,
                    "extensions": None,
                    "location_id": "ga4gh:VSL.wQ4TpNbsTPq_A-eQTL44gbP3f4fnp0vx",
                    "location": {
                        "id": None,
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 154170399},
                            "end": {"type": "Number", "value": 154170400},
                        },
                    },
                },
            },
            {
                "type": "GeneElement",
                "gene_descriptor": alk_gene_descriptor,
            },
            {
                "type": "LinkerSequenceElement",
                "linker_sequence": {
                    "id": "fusor.sequence:ACGT",
                    "type": "LiteralSequenceExpression",
                    "label": None,
                    "description": None,
                    "xrefs": None,
                    "alternate_labels": None,
                    "extensions": None,
                    "sequence_id": None,
                    "sequence": "ACGT",
                    "residue_type": "SO:0000348",
                },
            },
            {
                "type": "TemplatedSequenceElement",
                "region": {
                    "id": "fusor.location_descriptor:NC_000023.11",
                    "type": "LocationDescriptor",
                    "label": None,
                    "description": None,
                    "xrefs": None,
                    "alternate_labels": None,
                    "extensions": None,
                    "location_id": "ga4gh:VSL.zd12pX_ju2gLq9a9UOYgM8AtbkuhnyUu",
                    "location": {
                        "id": None,
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 44908820},
                            "end": {"type": "Number", "value": 44908822},
                        },
                    },
                },
                "strand": "+",
            },
            {"type": "MultiplePossibleGenesElement"},
        ],
        "regulatory_element": {
            "type": "RegulatoryElement",
            "regulatory_class": "promoter",
            "associatedGene": braf_gene_descriptor,
        },
    }


@pytest.fixture()
def fusion_example():
    """Create test fixture for a fake fusion without additional property expansion."""
    return {
        "type": "CategoricalFusion",
        "reading_frame_preserved": True,
        "critical_functional_domains": [
            {
                "type": "FunctionalDomain",
                "_id": "interpro:IPR020635",
                "label": "Tyrosine-protein kinase, catalytic domain",
                "status": "lost",
                "associatedGene": {
                    "id": "normalize.gene:hgnc%3A427",
                    "type": "Gene",
                    "label": "ALK",
                    "gene_id": "hgnc:427",
                },
                "sequence_location": {
                    "id": "fusor.location_descriptor:NP_004295.2",
                    "type": "LocationDescriptor",
                    "location": {
                        "type": "SequenceLocation",
                        "sequence_id": "refseq:NP_004295.2",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 1116},
                            "end": {"type": "Number", "value": 1383},
                        },
                    },
                },
            }
        ],
        "structure": [
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_152263.3",
                "exon_start": 1,
                "exon_start_offset": 0,
                "exon_end": 8,
                "exon_end_offset": 0,
                "gene_descriptor": {
                    "id": "normalize.gene:TPM3",
                    "type": "Gene",
                    "label": "TPM3",
                    "gene_id": "hgnc:12012",
                },
                "elementGenomicStart": {
                    "id": "fusor.location_descriptor:NC_000001.11",
                    "type": "LocationDescriptor",
                    "label": "NC_000001.11",
                    "location": {
                        "type": "SequenceLocation",
                        "sequence_id": "refseq:NC_000001.11",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 154192135},
                            "end": {"type": "Number", "value": 154192136},
                        },
                    },
                },
                "elementGenomicEnd": {
                    "id": "fusor.location_descriptor:NC_000001.11",
                    "type": "LocationDescriptor",
                    "label": "NC_000001.11",
                    "location": {
                        "type": "SequenceLocation",
                        "sequence_id": "refseq:NC_000001.11",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 154170399},
                            "end": {"type": "Number", "value": 154170400},
                        },
                    },
                },
            },
            {
                "type": "GeneElement",
                "gene_descriptor": {
                    "id": "normalize.gene:ALK",
                    "type": "Gene",
                    "label": "ALK",
                    "gene_id": "hgnc:427",
                },
            },
            {
                "type": "LinkerSequenceElement",
                "linker_sequence": {
                    "id": "fusor.sequence:ACGT",
                    "type": "LiteralSequenceExpression",
                    "sequence": "ACGT",
                    "residue_type": "SO:0000348",
                },
            },
            {
                "type": "TemplatedSequenceElement",
                "region": {
                    "id": "fusor.location_descriptor:NC_000023.11",
                    "type": "LocationDescriptor",
                    "location_id": "ga4gh:VSL.q0Hnb9gpYDyUuTix4Fesy5ungdnc4dWm",
                    "location": {
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 44908820},
                            "end": {"type": "Number", "value": 44908822},
                        },
                    },
                },
                "strand": "+",
            },
            {"type": "MultiplePossibleGenesElement"},
        ],
        "regulatory_element": {
            "type": "RegulatoryElement",
            "regulatory_class": "promoter",
            "associatedGene": {
                "id": "gene:BRAF",
                "type": "Gene",
                "label": "BRAF",
                "gene_id": "hgnc:1097",
            },
        },
    }
