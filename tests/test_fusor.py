"""Module for testing the FUSOR class."""
import pytest
from mock import patch
from fusor import FUSOR
from fusor.models import Fusion, GenomicRegionComponent, \
    TranscriptSegmentComponent


@pytest.fixture(scope="module")
def fusor():
    """Create test fixture for fusor object"""
    return FUSOR()


@pytest.fixture(scope="module")
def fusion():
    """Create fusion test fixture."""
    return {
        "r_frame_preserved": True,
        "protein_domains": [
            {
                "status": "lost",
                "name": "cystatin domain",
                "id": "interpro:IPR000010",
                "gene_descriptor": {
                    "id": "gene:CST1",
                    "gene_id": "hgnc:2743",
                    "label": "CST1",
                    "type": "GeneDescriptor"
                }
            }
        ],
        "structural_components": [
            {
                "component_type": "transcript_segment",
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
                "component_genomic_region": {
                    "id": "refseq:NM_152263.3_exon1-exon8",
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
                                "value": 154170399
                            },
                            "type": "SequenceInterval"
                        }
                    }
                }
            },
            {
                "component_type": "gene",
                "gene_descriptor": {
                    "id": "gene:ALK",
                    "type": "GeneDescriptor",
                    "gene_id": "hgnc:427",
                    "label": "ALK"
                }
            },
            {
                "component_type": "linker_sequence",
                "linker_sequence": {
                    "id": "sequence:ACGT",
                    "type": "SequenceDescriptor",
                    "sequence": "ACGT",
                    "residue_type": "SO:0000348"
                }
            },
            {
                "component_type": "genomic_region",
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
                "component_type": "unknown_gene"
            }
        ],
        "causative_event": "rearrangement",
        "regulatory_elements": [
            {
                "type": "promoter",
                "gene_descriptor": {
                    "id": "gene:BRAF",
                    "type": "GeneDescriptor",
                    "gene_id": "hgnc:1097",
                    "label": "BRAF"
                }
            }
        ]
    }


SEQUENCE_IDS = ["ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT",
                "ga4gh:SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl"]


@patch.object(FUSOR, "translate_identifier",
              side_effect=SEQUENCE_IDS)
def test_add_additional_fields(test_translate_identifier, fusor,
                               fusion_example, fusion):
    """Test that add_additional_fields method works correctly."""
    expected_fusion = Fusion(**fusion_example)
    fusion = Fusion(**fusion)
    fusor.add_additional_fields(fusion)
    assert fusion.dict() == expected_fusion.dict()


@patch.object(FUSOR, "translate_identifier",
              side_effect=SEQUENCE_IDS)
def test_add_sequence_id(test_translate_identifier, fusor, fusion_example,
                         fusion):
    """Test that add_sequence_id method works correctly."""
    expected_fusion = Fusion(**fusion_example)
    test_fusion = Fusion(**fusion)
    for structural_component in expected_fusion.structural_components:
        if isinstance(structural_component, GenomicRegionComponent):
            structural_component.region.location_id = None
        elif isinstance(structural_component, TranscriptSegmentComponent):
            structural_component.component_genomic_region.location_id = None
    fusor.add_sequence_id(test_fusion)
    assert test_fusion.dict() == expected_fusion.dict()


def test_add_location_id(fusor, fusion_example, fusion):
    """Test that add_location_id method works correctly."""
    expected_fusion = Fusion(**fusion_example)
    test_fusion = Fusion(**fusion)
    fusor.add_sequence_id(test_fusion)
    fusor.add_location_id(test_fusion)
    assert test_fusion.dict() == expected_fusion.dict()
