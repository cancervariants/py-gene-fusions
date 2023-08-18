"""Module containing methods and fixtures used throughout tests."""
import pytest

from fusor.fusor import FUSOR


@pytest.fixture(scope="session")
def fusor_instance():
    """Create test fixture for fusor object"""
    return FUSOR()


@pytest.fixture(scope="session")
def braf_gene_descriptor():
    """Create gene descriptor params for BRAF."""
    return {
        "id": "normalize.gene:BRAF",
        "type": "GeneDescriptor",
        "label": "BRAF",
        "xrefs": ["ensembl:ENSG00000157764", "ncbigene:673"],
        "alternate_labels": ["BRAF1", "BRAF-1", "NS7", "B-raf", "B-RAF1", "RAFB1"],
        "extensions": [
            {"type": "Extension", "name": "symbol_status", "value": "approved"},
            {
                "type": "Extension",
                "name": "approved_name",
                "value": "B-Raf proto-oncogene, serine/threonine kinase",
            },
            {
                "type": "Extension",
                "name": "hgnc_locations",
                "value": [
                    {
                        "_id": "ga4gh:VCL.O6yCQ1cnThOrTfK9YUgMlTfM6HTqbrKw",
                        "type": "ChromosomeLocation",
                        "species_id": "taxonomy:9606",
                        "chr": "7",
                        "interval": {
                            "end": "q34",
                            "start": "q34",
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
                        "_id": "ga4gh:VSL.amNWL6i7F2nbSZAf2QLTRTujxuDrd0pR",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        "interval": {
                            "start": {"type": "Number", "value": 140719326},
                            "end": {"type": "Number", "value": 140924929},
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
                        "_id": "ga4gh:VCL.O6yCQ1cnThOrTfK9YUgMlTfM6HTqbrKw",
                        "type": "ChromosomeLocation",
                        "species_id": "taxonomy:9606",
                        "chr": "7",
                        "interval": {
                            "end": "q34",
                            "start": "q34",
                            "type": "CytobandInterval",
                        },
                    },
                    {
                        "_id": "ga4gh:VSL.xZU3kL8F6t2ca6WH_26CWKfNW9-owhR4",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        "interval": {
                            "start": {"type": "Number", "value": 140713327},
                            "end": {"type": "Number", "value": 140924929},
                            "type": "SequenceInterval",
                        },
                    },
                ],
            },
            {
                "type": "Extension",
                "name": "associated_with",
                "value": [
                    "pubmed:2284096",
                    "refseq:NM_004333",
                    "iuphar:1943",
                    "orphanet:119066",
                    "cosmic:BRAF",
                    "ena.embl:M95712",
                    "ccds:CCDS87555",
                    "ucsc:uc003vwc.5",
                    "pubmed:1565476",
                    "vega:OTTHUMG00000157457",
                    "uniprot:P15056",
                    "ccds:CCDS5863",
                    "ccds:CCDS94218",
                    "ccds:CCDS94219",
                    "omim:164757",
                ],
            },
            {
                "type": "Extension",
                "name": "hgnc_locus_type",
                "value": "gene with protein product",
            },
            {"type": "Extension", "name": "ncbi_gene_type", "value": "protein-coding"},
            {"type": "Extension", "name": "ensembl_biotype", "value": "protein_coding"},
        ],
        "gene_id": "hgnc:1097",
    }


@pytest.fixture(scope="session")
def alk_gene_descriptor():
    """Create test fixture for ALK gene descriptor params"""
    return {
        "id": "normalize.gene:ALK",
        "type": "GeneDescriptor",
        "label": "ALK",
        "description": None,
        "xrefs": ["ensembl:ENSG00000171094", "ncbigene:238"],
        "alternate_labels": ["NBLST3", "CD246", "ALK1"],
        "extensions": [
            {
                "type": "Extension",
                "name": "symbol_status",
                "value": "approved",
            },
            {
                "type": "Extension",
                "name": "approved_name",
                "value": "ALK receptor tyrosine kinase",
            },
            {
                "type": "Extension",
                "name": "hgnc_locations",
                "value": [
                    {
                        "_id": "ga4gh:VCL.VE7uJHat7zIWFf9AzNM85jj05r1dLzsD",
                        "type": "ChromosomeLocation",
                        "species_id": "taxonomy:9606",
                        "chr": "2",
                        "interval": {
                            "end": "p23.1",
                            "start": "p23.2",
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
                        "_id": "ga4gh:VSL.-k3kxW3qMyV-oBTvTffVZojkJBLs0flu",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
                        "interval": {
                            "start": {"type": "Number", "value": 29192773},
                            "end": {"type": "Number", "value": 29921586},
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
                        "_id": "ga4gh:VCL.VE7uJHat7zIWFf9AzNM85jj05r1dLzsD",
                        "type": "ChromosomeLocation",
                        "species_id": "taxonomy:9606",
                        "chr": "2",
                        "interval": {
                            "end": "p23.1",
                            "start": "p23.2",
                            "type": "CytobandInterval",
                        },
                    },
                    {
                        "_id": "ga4gh:VSL.-k3kxW3qMyV-oBTvTffVZojkJBLs0flu",
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.pnAqCRBrTsUoBghSD1yp_jXWSmlbdh4g",
                        "interval": {
                            "start": {"type": "Number", "value": 29192773},
                            "end": {"type": "Number", "value": 29921586},
                            "type": "SequenceInterval",
                        },
                    },
                ],
            },
            {
                "type": "Extension",
                "name": "associated_with",
                "value": [
                    "ccds:CCDS33172",
                    "pubmed:8122112",
                    "orphanet:160020",
                    "ccds:CCDS86828",
                    "cosmic:ALK",
                    "uniprot:Q9UM73",
                    "omim:105590",
                    "iuphar:1839",
                    "hcdmdb:CD246",
                    "vega:OTTHUMG00000152034",
                    "ena.embl:D45915",
                    "refseq:NM_004304",
                    "ucsc:uc002rmy.4",
                ],
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
        ],
        "gene_id": "hgnc:427",
        "gene": None,
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
                "associated_gene": alk_gene_descriptor,
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
        "structural_elements": [
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_152263.3",
                "exon_start": 1,
                "exon_start_offset": 0,
                "exon_end": 8,
                "exon_end_offset": 0,
                "gene_descriptor": {
                    "id": "normalize.gene:TPM3",
                    "type": "GeneDescriptor",
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
                        "CMYP4A",
                        "CMYP4B",
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
                                    "sequence_id": "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",  # noqa: E501
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
                                    "sequence_id": "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",  # noqa: E501
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
                                "ccds:CCDS91060",
                                "ccds:CCDS91061",
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
                    ],
                    "gene_id": "hgnc:12012",
                    "gene": None,
                },
                "element_genomic_start": {
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
                "element_genomic_end": {
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
                    "type": "SequenceDescriptor",
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
            "associated_gene": braf_gene_descriptor,
        },
    }


@pytest.fixture(scope="function")
def fusion_example():
    """Create test fixture for a fake fusion without additional property expansion."""
    return {
        "type": "CategoricalFusion",
        "r_frame_preserved": True,
        "critical_functional_domains": [
            {
                "type": "FunctionalDomain",
                "_id": "interpro:IPR020635",
                "label": "Tyrosine-protein kinase, catalytic domain",
                "status": "lost",
                "associated_gene": {
                    "id": "normalize.gene:hgnc%3A427",
                    "type": "GeneDescriptor",
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
        "structural_elements": [
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_152263.3",
                "exon_start": 1,
                "exon_start_offset": 0,
                "exon_end": 8,
                "exon_end_offset": 0,
                "gene_descriptor": {
                    "id": "normalize.gene:TPM3",
                    "type": "GeneDescriptor",
                    "label": "TPM3",
                    "gene_id": "hgnc:12012",
                },
                "element_genomic_start": {
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
                "element_genomic_end": {
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
                    "type": "GeneDescriptor",
                    "label": "ALK",
                    "gene_id": "hgnc:427",
                },
            },
            {
                "type": "LinkerSequenceElement",
                "linker_sequence": {
                    "id": "fusor.sequence:ACGT",
                    "type": "SequenceDescriptor",
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
            "associated_gene": {
                "id": "gene:BRAF",
                "type": "GeneDescriptor",
                "label": "BRAF",
                "gene_id": "hgnc:1097",
            },
        },
    }
