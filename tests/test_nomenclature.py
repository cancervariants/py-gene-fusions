"""Test nomenclature generation."""
import pytest

from fusor.models import CategoricalFusion, AssayedFusion


@pytest.fixture(scope="module")
def reg_example():
    """Nonsense fusion testing correct regulatory element description."""
    return AssayedFusion(
        **{
            "type": "AssayedFusion",
            "regulatory_element": {
                "type": "RegulatoryElement",
                "regulatory_class": "riboswitch",
                "associated_gene": {
                    "id": "normalize.gene:ABL1",
                    "type": "GeneDescriptor",
                    "label": "ABL1",
                    "gene_id": "hgnc:76",
                },
            },
            "structural_elements": [
                {
                    "type": "GeneElement",
                    "gene_descriptor": {
                        "id": "normalize.gene:BCR",
                        "type": "GeneDescriptor",
                        "label": "BCR",
                        "gene_id": "hgnc:1014",
                    },
                },
                {"type": "UnknownGeneElement"},
            ],
            "causative_event": {
                "type": "CausativeEvent",
                "event_type": "rearrangement",
            },
            "assay": {
                "type": "Assay",
                "assay_name": "a",
                "assay_id": "a:b",
                "method_uri": "a:b",
                "fusion_detection": "observed",
            },
        }
    )


@pytest.fixture(scope="module")
def reg_location_example():
    """Nonsense fusion testing correct regulatory element description."""
    return AssayedFusion(
        **{
            "type": "AssayedFusion",
            "regulatory_element": {
                "type": "RegulatoryElement",
                "regulatory_class": "promoter",
                "associated_gene": {
                    "id": "normalize.gene:P2RY8",
                    "type": "GeneDescriptor",
                    "label": "P2RY8",
                    "gene_id": "hgnc:15524",
                },
                "feature_location": {
                    "type": "LocationDescriptor",
                    "id": "fusor.location_descriptor:NC_000023.11",
                    "location": {
                        "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                        "type": "SequenceLocation",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 1462581},
                            "end": {"type": "Number", "value": 1534182},
                        },
                    },
                },
            },
            "structural_elements": [
                {
                    "type": "GeneElement",
                    "gene_descriptor": {
                        "id": "normalize.gene:SOX5",
                        "type": "GeneDescriptor",
                        "label": "SOX5",
                        "gene_id": "hgnc:11201",
                    },
                },
            ],
            "causative_event": {
                "type": "CausativeEvent",
                "event_type": "rearrangement",
            },
            "assay": {
                "type": "Assay",
                "assay_name": "a",
                "assay_id": "a:b",
                "method_uri": "a:b",
                "fusion_detection": "observed",
            },
        }
    )


@pytest.fixture(scope="module")
def exon_offset_example():
    """Provide example of tx segment with positive exon end offset"""
    return CategoricalFusion(
        **{
            "type": "CategoricalFusion",
            "structural_elements": [
                {
                    "type": "GeneElement",
                    "gene_descriptor": {
                        "id": "normalize.gene:BRAF",
                        "type": "GeneDescriptor",
                        "label": "BRAF",
                        "gene_id": "hgnc:1097",
                    },
                },
                {
                    "type": "TranscriptSegmentElement",
                    "transcript": "refseq:NM_002529.3",
                    "exon_start": 2,
                    "exon_start_offset": 20,
                    "gene_descriptor": {
                        "id": "normalize.gene:NTRK1",
                        "type": "GeneDescriptor",
                        "label": "NTRK1",
                        "gene_id": "hgnc:8031",
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
                                "start": {"type": "Number", "value": 156864428},
                                "end": {"type": "Number", "value": 156864429},
                            },
                        },
                    },
                },
            ],
        }
    )


def test_generate_nomenclature(
    fusor_instance,
    fusion_example,
    exhaustive_example,
    reg_example,
    reg_location_example,
    exon_offset_example,
):
    """Test that nomenclature generation is correct."""
    fixture_nomenclature = "reg_p@BRAF(hgnc:1097)::refseq:NM_152263.3(TPM3):e.1_8::ALK(hgnc:427)::ACGT::refseq:NC_000023.11(chr X):g.44908820_44908822(+)::v"  # noqa: E501
    nm = fusor_instance.generate_nomenclature(CategoricalFusion(**fusion_example))
    assert nm == fixture_nomenclature
    nm = fusor_instance.generate_nomenclature(CategoricalFusion(**exhaustive_example))
    assert nm == fixture_nomenclature

    from fusor import examples

    nm = fusor_instance.generate_nomenclature(examples.bcr_abl1)
    assert (nm == "refseq:NM_004327.3(BCR):e.2+182::ACTAAAGCG::refseq:NM_005157.5(ABL1):e.2-173")  # noqa: E501

    nm = fusor_instance.generate_nomenclature(examples.bcr_abl1_expanded)
    assert (nm == "refseq:NM_004327.3(BCR):e.2+182::ACTAAAGCG::refseq:NM_005157.5(ABL1):e.2-173")  # noqa: E501

    nm = fusor_instance.generate_nomenclature(examples.alk)
    assert nm == "ALK(hgnc:427)::v"

    nm = fusor_instance.generate_nomenclature(examples.tpm3_ntrk1)
    assert (
        nm == "refseq:NM_152263.3(TPM3):e.8::refseq:NM_002529.3(NTRK1):e.10"
    )  # noqa: E501

    nm = fusor_instance.generate_nomenclature(examples.tpm3_pdgfrb)
    assert (
        nm == "refseq:NM_152263.3(TPM3):e.8::refseq:NM_002609.3(PDGFRB):e.11"
    )  # noqa: E501

    nm = fusor_instance.generate_nomenclature(examples.ewsr1)
    assert nm == "EWSR1(hgnc:3508)::?"

    nm = fusor_instance.generate_nomenclature(examples.igh_myc)
    assert nm == "reg_e_EH38E3121735@IGH(hgnc:5477)::MYC(hgnc:7553)"

    nm = fusor_instance.generate_nomenclature(reg_example)
    assert nm == "reg_riboswitch@ABL1(hgnc:76)::BCR(hgnc:1014)::?"

    nm = fusor_instance.generate_nomenclature(reg_location_example)
    assert (
        nm
        == "reg_p_refseq:NC_000023.11(chr X):g.1462581_1534182@P2RY8(hgnc:15524)::SOX5(hgnc:11201)"  # noqa: E501
    )

    nm = fusor_instance.generate_nomenclature(exon_offset_example)
    assert nm == "BRAF(hgnc:1097)::refseq:NM_002529.3(NTRK1):e.2+20"
