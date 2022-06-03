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


@pytest.fixture(scope="module")
def fusion():
    """Create fusion test fixture."""
    return {
        "type": "CategoricalFusion",
        "r_frame_preserved": True,
        "critical_functional_domains": [
            {
                "status": "lost",
                "label": "Tyrosine-protein kinase, catalytic domain",
                "id": "interpro:IPR020635",
                "associated_gene": {
                    "id": "gene:ALK",
                    "gene_id": "hgnc:1837",
                    "label": "ALK"
                },
                "sequence_location": {
                    "id": "fusor.location_descriptor:NP_002520.2",
                    "type": "LocationDescriptor",
                    "location": {
                        "sequence_id": "ga4gh:SQ.vJvm06Wl5J7DXHynR9ksW7IK3_3jlFK6",  # noqa: E501
                        "type": "SequenceLocation",
                        "interval": {
                            "start": {
                                "type": "Number",
                                "value": 510
                            },
                            "end": {
                                "type": "Number",
                                "value": 781
                            }
                        }
                    }
                }
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
                    "id": "gene:TPM3",
                    "gene_id": "hgnc:12012",
                    "type": "GeneDescriptor",
                    "label": "TPM3"
                },
                "element_genomic_start": {
                    "id": "refseq:NM_152263.3_exon1",
                    "type": "LocationDescriptor",
                    "location": {
                        "sequence_id": "refseq:NM_152263.3",
                        "type": "SequenceLocation",
                        "interval": {
                            "start": {
                                "type": "Number",
                                "value": 154192135
                            },
                            "end": {
                                "type": "Number",
                                "value": 154192136
                            },
                            "type": "SequenceInterval"
                        }
                    }
                },
                "element_genomic_end": {
                    "id": "refseq:NM_152263.3_exon8",
                    "type": "LocationDescriptor",
                    "location": {
                        "sequence_id": "refseq:NM_152263.3",
                        "type": "SequenceLocation",
                        "interval": {
                            "start": {
                                "type": "Number",
                                "value": 154170398
                            },
                            "end": {
                                "type": "Number",
                                "value": 154170399
                            },
                            "type": "SequenceInterval"
                        }
                    }
                }
            },
            {
                "type": "GeneElement",
                "gene_descriptor": {
                    "id": "gene:ALK",
                    "type": "GeneDescriptor",
                    "gene_id": "hgnc:427",
                    "label": "ALK"
                }
            },
            {
                "type": "LinkerSequenceElement",
                "linker_sequence": {
                    "id": "sequence:ACGT",
                    "type": "SequenceDescriptor",
                    "sequence": "ACGT",
                    "residue_type": "SO:0000348"
                }
            },
            {
                "type": "TemplatedSequenceElement",
                "region": {
                    "id": "chr12:44908821-44908822(+)",
                    "type": "LocationDescriptor",
                    "location": {
                        "type": "SequenceLocation",
                        "sequence_id": "refseq:NC_000012.12",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {
                                "type": "Number",
                                "value": 44908821
                            },
                            "end": {
                                "type": "Number",
                                "value": 44908822
                            }
                        }
                    },
                    "label": "chr12:44908821-44908822(+)"
                },
                "strand": "+"
            },
            {
                "type": "MultiplePossibleGenesElement"
            }
        ],
        "regulatory_elements": [
            {
                "type": "RegulatoryElement",
                "regulatory_class": "promoter",
                "associated_gene": {
                    "id": "gene:BRAF",
                    "type": "GeneDescriptor",
                    "gene_id": "hgnc:1097",
                    "label": "BRAF"
                }
            }
        ]
    }
