"""Test nomenclature generation."""

import pytest
from ga4gh.core.domain_models import Gene

from fusor.models import AssayedFusion, CategoricalFusion, TranscriptSegmentElement
from fusor.nomenclature import tx_segment_nomenclature


@pytest.fixture(scope="module")
def reg_example():
    """Nonsense fusion testing correct regulatory element description."""
    return AssayedFusion(
        type="AssayedFusion",
        regulatoryElement={
            "type": "RegulatoryElement",
            "regulatoryClass": "riboswitch",
            "associatedGene": {
                "type": "Gene",
                "label": "ABL1",
                "id": "hgnc:76",
            },
        },
        structure=[
            {
                "type": "GeneElement",
                "gene": {
                    "type": "Gene",
                    "label": "BCR",
                    "id": "hgnc:1014",
                },
            },
            {"type": "UnknownGeneElement"},
        ],
        causativeEvent={
            "type": "CausativeEvent",
            "eventType": "rearrangement",
        },
        assay={
            "type": "Assay",
            "assayName": "a",
            "assayId": "a:b",
            "methodUti": "a:b",
            "fusionDetection": "observed",
        },
    )


@pytest.fixture(scope="module")
def reg_location_example():
    """Nonsense fusion testing correct regulatory element description."""
    return AssayedFusion(
        type="AssayedFusion",
        regulatoryElement={
            "type": "RegulatoryElement",
            "regulatoryClass": "promoter",
            "associatedGene": {
                "type": "Gene",
                "label": "P2RY8",
                "id": "hgnc:15524",
            },
            "featureLocation": {
                "type": "SequenceLocation",
                "id": "ga4gh:SL.KMHXvX8m5fD8PcGlQu2Vja3m7bt2iqfK",
                "sequenceReference": {
                    "id": "refseq:NC_000023.11",
                    "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                },
                "start": 1462581,
                "end": 1534182,
            },
        },
        structural_elements=[
            {
                "type": "GeneElement",
                "gene": {
                    "type": "Gene",
                    "label": "SOX5",
                    "id": "hgnc:11201",
                },
            },
        ],
        causative_event={
            "type": "CausativeEvent",
            "eventType": "rearrangement",
        },
        assay={
            "type": "Assay",
            "assayName": "a",
            "assayId": "a:b",
            "methodUri": "a:b",
            "fusionDetection": "observed",
        },
    )


@pytest.fixture(scope="module")
def exon_offset_example():
    """Provide example of tx segment with positive exon end offset"""
    return CategoricalFusion(
        type="CategoricalFusion",
        structural_elements=[
            {
                "type": "GeneElement",
                "gene": {
                    "type": "Gene",
                    "label": "BRAF",
                    "id": "hgnc:1097",
                },
            },
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_002529.3",
                "exonStart": 2,
                "exonStartOffset": 20,
                "gene": {
                    "type": "Gene",
                    "label": "NTRK1",
                    "id": "hgnc:8031",
                },
                "elementGenomicStart": {
                    "id": "ga4gh:SL.XEvDpRaKgoeQuQrhRwGzGK2uanHY4en8",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "refseq:NC_000001.11",
                        "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                    },
                    "start": {"type": "Number", "value": 156864428},
                    "end": {"type": "Number", "value": 156864429},
                },
            },
        ],
    )


@pytest.fixture(scope="module")
def tx_seg_example():
    """Provide example of transcript segment element."""
    return TranscriptSegmentElement(
        type="TranscriptSegmentElement",
        transcript="refseq:NM_152263.3",
        exon_start=1,
        exon_start_offset=0,
        exon_end=8,
        exon_end_offset=0,
        gene=Gene(id="hgnc:12012", label="TPM3"),
        element_genomic_start={
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
        element_genomic_end={
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
    )


@pytest.fixture(scope="module")
def junction_example():
    """Provide example of junction element."""
    return TranscriptSegmentElement(
        type="TranscriptSegmentElement",
        transcript="refseq:NM_152263.3",
        exon_end=8,
        exon_end_offset=0,
        gene_descriptor={
            "id": "normalize.gene:TPM3",
            "type": "GeneDescriptor",
            "label": "TPM3",
            "gene_id": "hgnc:12012",
        },
        element_genomic_end={
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
    fixture_nomenclature = "reg_p@BRAF(hgnc:1097)::NM_152263.3(TPM3):e.1_8::ALK(hgnc:427)::ACGT::NC_000023.11(chr X):g.44908820_44908822(+)::v"
    nm = fusor_instance.generate_nomenclature(CategoricalFusion(**fusion_example))
    assert nm == fixture_nomenclature
    nm = fusor_instance.generate_nomenclature(CategoricalFusion(**exhaustive_example))
    assert nm == fixture_nomenclature

    from fusor import examples

    nm = fusor_instance.generate_nomenclature(examples.bcr_abl1)
    assert nm == "NM_004327.3(BCR):e.2+182::ACTAAAGCG::NM_005157.5(ABL1):e.2-173"

    nm = fusor_instance.generate_nomenclature(examples.bcr_abl1_expanded)
    assert nm == "NM_004327.3(BCR):e.2+182::ACTAAAGCG::NM_005157.5(ABL1):e.2-173"

    nm = fusor_instance.generate_nomenclature(examples.alk)
    assert nm == "v::ALK(hgnc:427)"

    nm = fusor_instance.generate_nomenclature(examples.tpm3_ntrk1)
    assert nm == "NM_152263.3(TPM3):e.8(::)NM_002529.3(NTRK1):e.10"

    nm = fusor_instance.generate_nomenclature(examples.tpm3_pdgfrb)
    assert nm == "NM_152263.3(TPM3):e.1_8::NM_002609.3(PDGFRB):e.11_22"

    nm = fusor_instance.generate_nomenclature(examples.ewsr1)
    assert nm == "EWSR1(hgnc:3508)(::)?"

    nm = fusor_instance.generate_nomenclature(examples.ewsr1_no_assay)
    assert nm == "EWSR1(hgnc:3508)::?"

    nm = fusor_instance.generate_nomenclature(examples.ewsr1_no_causative_event)
    assert nm == "EWSR1(hgnc:3508)(::)?"

    nm = fusor_instance.generate_nomenclature(examples.ewsr1_elements_only)
    assert nm == "EWSR1(hgnc:3508)::?"

    nm = fusor_instance.generate_nomenclature(examples.igh_myc)
    assert nm == "reg_e_EH38E3121735@IGH(hgnc:5477)::MYC(hgnc:7553)"

    nm = fusor_instance.generate_nomenclature(reg_example)
    assert nm == "reg_riboswitch@ABL1(hgnc:76)::BCR(hgnc:1014)::?"

    nm = fusor_instance.generate_nomenclature(reg_location_example)
    assert (
        nm
        == "reg_p_NC_000023.11(chr X):g.1462581_1534182@P2RY8(hgnc:15524)::SOX5(hgnc:11201)"
    )

    nm = fusor_instance.generate_nomenclature(exon_offset_example)
    assert nm == "BRAF(hgnc:1097)::NM_002529.3(NTRK1):e.2+20"


def test_component_nomenclature(tx_seg_example, junction_example):
    """Test that individual object nomenclature generators are correct."""
    nm = tx_segment_nomenclature(tx_seg_example)
    assert nm == "NM_152263.3(TPM3):e.1_8"

    nm = tx_segment_nomenclature(junction_example)
    assert nm == "NM_152263.3(TPM3):e.8"
