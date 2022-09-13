"""Module containing methods and fixtures used throughout tests."""
import pytest

from fusor.fusor import FUSOR


@pytest.fixture(scope="session")
def fusor_instance():
    """Create test fixture for fusor object"""
    return FUSOR()


@pytest.fixture(scope="module")
def exhaustive_example():
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
                "associated_gene": {
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
                            "name": "chromosome_location",
                            "value": {
                                "_id": "ga4gh:VCL.VE7uJHat7zIWFf9AzNM85jj05r1dLzsD",
                                "type": "ChromosomeLocation",
                                "species_id": "taxonomy:9606",
                                "chr": "2",
                                "interval": {
                                    "type": "CytobandInterval",
                                    "start": "p23.2",
                                    "end": "p23.1",
                                },
                            },
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
                },
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
                            "name": "chromosome_location",
                            "value": {
                                "_id": "ga4gh:VCL.rmJvYV5JccRSEoMVxe5BmuHs9S2VZ4uR",
                                "type": "ChromosomeLocation",
                                "species_id": "taxonomy:9606",
                                "chr": "1",
                                "interval": {
                                    "type": "CytobandInterval",
                                    "start": "q21.3",
                                    "end": "q21.3",
                                },
                            },
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
                "gene_descriptor": {
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
                            "name": "chromosome_location",
                            "value": {
                                "_id": "ga4gh:VCL.VE7uJHat7zIWFf9AzNM85jj05r1dLzsD",
                                "type": "ChromosomeLocation",
                                "species_id": "taxonomy:9606",
                                "chr": "2",
                                "interval": {
                                    "type": "CytobandInterval",
                                    "start": "p23.2",
                                    "end": "p23.1",
                                },
                            },
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
                },
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
            "associated_gene": {
                "id": "normalize.gene:BRAF",
                "type": "GeneDescriptor",
                "label": "BRAF",
                "description": None,
                "xrefs": ["ensembl:ENSG00000157764", "ncbigene:673"],
                "alternate_labels": [
                    "BRAF1",
                    "BRAF-1",
                    "NS7",
                    "B-raf",
                    "B-RAF1",
                    "RAFB1",
                ],
                "extensions": [
                    {"type": "Extension", "name": "symbol_status", "value": "approved"},
                    {
                        "type": "Extension",
                        "name": "approved_name",
                        "value": "B-Raf proto-oncogene, serine/threonine kinase",
                    },
                    {
                        "type": "Extension",
                        "name": "chromosome_location",
                        "value": {
                            "_id": "ga4gh:VCL.O6yCQ1cnThOrTfK9YUgMlTfM6HTqbrKw",
                            "type": "ChromosomeLocation",
                            "species_id": "taxonomy:9606",
                            "chr": "7",
                            "interval": {
                                "type": "CytobandInterval",
                                "start": "q34",
                                "end": "q34",
                            },
                        },
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
                            "omim:164757",
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
                "gene_id": "hgnc:1097",
                "gene": None,
            },
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
