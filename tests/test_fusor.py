"""Module for testing the FUSOR class."""
import pytest
from ga4gh.vrsatile.pydantic.vrsatile_model import GeneDescriptor
from fusor import FUSOR
from fusor.models import Fusion, TemplatedSequenceComponent, \
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
                    "gene_id": "hgnc:2473",
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
                "component_genomic_start": {
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
                "component_genomic_end": {
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
                "component_type": "templated_sequence",
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


def compare_gene_descriptor(actual, expected):
    """Test that actual and expected gene descriptors match."""
    assert actual["id"] == expected["id"]
    assert actual["type"] == expected["type"]
    assert actual["gene_id"] == expected["gene_id"]
    assert actual["label"] == expected["label"]
    assert set(actual["xrefs"]) == set(expected["xrefs"]), "xrefs"
    if expected["alternate_labels"]:
        assert set(actual["alternate_labels"]) == set(expected["alternate_labels"]), "alt labels"  # noqa: E501
    else:
        assert actual["alternate_labels"] == expected["alternate_labels"]
    extensions_present = "extensions" in expected.keys()
    assert ("extensions" in actual.keys()) == extensions_present
    if extensions_present:
        assert len(actual["extensions"]) == len(expected["extensions"]), \
            "len of extensions"
        n_ext_correct = 0
        for expected_ext in expected["extensions"]:
            for actual_ext in actual["extensions"]:
                if actual_ext["name"] == expected_ext["name"]:
                    assert isinstance(actual_ext["value"],
                                      type(expected_ext["value"]))
                    if isinstance(expected_ext["value"], list):
                        assert set(actual_ext["value"]) == \
                               set(expected_ext["value"]), f"{expected_ext['value']} value"  # noqa: E501
                    else:
                        assert actual_ext["value"] == expected_ext["value"]
                    assert actual_ext["type"] == expected_ext["type"]
                    n_ext_correct += 1
        assert n_ext_correct == len(expected["extensions"]), \
            "number of correct extensions"


def test_add_additional_fields(fusor, fusion_example, fusion):
    """Test that add_additional_fields method works correctly."""
    expected_fusion = Fusion(**fusion_example)
    fusion = Fusion(**fusion)
    fusor.add_additional_fields(fusion)
    assert fusion.dict() == expected_fusion.dict()


def test_add_sequence_id(fusor, fusion_example, fusion):
    """Test that add_sequence_id method works correctly."""
    expected_fusion = Fusion(**fusion_example)
    test_fusion = Fusion(**fusion)
    for structural_component in expected_fusion.structural_components:
        if isinstance(structural_component, TemplatedSequenceComponent):
            structural_component.region.location_id = None
        elif isinstance(structural_component, TranscriptSegmentComponent):
            if structural_component.component_genomic_start:
                structural_component.component_genomic_start.location_id = None
            if structural_component.component_genomic_end:
                structural_component.component_genomic_end.location_id = None
    fusor.add_sequence_id(test_fusion)
    assert test_fusion.dict() == expected_fusion.dict()


def test_add_location_id(fusor, fusion_example, fusion):
    """Test that add_location_id method works correctly."""
    expected_fusion = Fusion(**fusion_example)
    test_fusion = Fusion(**fusion)
    fusor.add_sequence_id(test_fusion)
    fusor.add_location_id(test_fusion)
    assert test_fusion.dict() == expected_fusion.dict()


def test_translate_identifier(fusor):
    """Test that translate_identifier method works correctly."""
    expected = "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT"
    identifier = fusor.translate_identifier("NM_152263.3")
    assert identifier == expected

    identifier = fusor.translate_identifier("refseq:NM_152263.3")
    assert identifier == expected

    identifier = fusor.translate_identifier("refseq_152263.3")
    assert identifier is None


def test__normalized_gene_descriptor(fusor):
    """Test that _normalized_gene_descriptor works correctly."""
    # Actual response is tested in test_add_gene_descriptor
    resp = fusor._normalized_gene_descriptor("BRAF")
    assert resp[0]
    assert resp[1] is None
    assert isinstance(resp[0], GeneDescriptor)

    resp = fusor._normalized_gene_descriptor("B R A F")
    assert resp[0] is None
    assert resp[1] == "gene-normalizer unable to normalize B R A F"


def test_add_gene_descriptor(fusor, exhaustive_example, fusion):
    """Test that add_gene_descriptor method works correctly."""
    expected_fusion = Fusion(**exhaustive_example)
    test_fusion = Fusion(**fusion)
    fusor.add_sequence_id(test_fusion)
    fusor.add_location_id(test_fusion)
    fusor.add_gene_descriptor(test_fusion)

    e_gds = set()
    t_gds = set()
    for e_field in [expected_fusion.protein_domains,
                    expected_fusion.structural_components,
                    expected_fusion.regulatory_elements]:
        for t_field in [test_fusion.protein_domains,
                        test_fusion.structural_components,
                        test_fusion.regulatory_elements]:
            for e_obj in e_field:
                for t_obj in t_field:
                    if "gene_descriptor" in e_obj.__fields__.keys():
                        e_gd = e_obj.gene_descriptor.label
                        e_gds.add(e_gd)
                        if "gene_descriptor" in t_obj.__fields__.keys():
                            t_gd = t_obj.gene_descriptor.label
                            t_gds.add(t_gd)
                            if e_gd == t_gd:
                                compare_gene_descriptor(
                                    t_obj.gene_descriptor.dict(),
                                    e_obj.gene_descriptor.dict()
                                )
    assert t_gds == e_gds
