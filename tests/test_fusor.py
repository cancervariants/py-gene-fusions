"""Module for testing the FUSOR class."""
import copy

import pytest
from ga4gh.vrsatile.pydantic.vrsatile_models import GeneDescriptor, LocationDescriptor

from fusor.exceptions import FUSORParametersException
from fusor.models import (
    AssayedFusion,
    CategoricalFusion,
    FunctionalDomain,
    GeneElement,
    LinkerElement,
    MultiplePossibleGenesElement,
    RegulatoryClass,
    RegulatoryElement,
    TemplatedSequenceElement,
    TranscriptSegmentElement,
    UnknownGeneElement,
)


@pytest.fixture(scope="module")
def braf_gene_descr_min():
    """Create minimal gene descriptor for BRAF"""
    return GeneDescriptor(id="normalize.gene:BRAF", label="BRAF", gene_id="hgnc:1097")


@pytest.fixture(scope="module")
def braf_gene_descr(braf_gene_descriptor):
    """Create gene descriptor object for braf"""
    return GeneDescriptor(**braf_gene_descriptor)


@pytest.fixture(scope="module")
def linker_element():
    """Create linker element test fixture."""
    params = {
        "linker_sequence": {
            "id": "fusor.sequence:ACT",
            "sequence": "ACT",
            "residue_type": "SO:0000348",
            "type": "SequenceDescriptor",
        },
        "type": "LinkerSequenceElement",
    }
    return LinkerElement(**params)


@pytest.fixture(scope="module")
def location_descriptor_braf_domain():
    """Create location descriptor fixture for BRAF catalytic domain"""
    params = {
        "id": "fusor.location_descriptor:NP_004324.2",
        "type": "LocationDescriptor",
        "location": {
            "sequence_id": "refseq:NP_004324.2",
            "type": "SequenceLocation",
            "interval": {
                "start": {"type": "Number", "value": 458},
                "end": {
                    "type": "Number",
                    "value": 712,
                },
            },
        },
    }
    return LocationDescriptor(**params)


@pytest.fixture(scope="module")
def location_descriptor_braf_domain_seq_id():
    """Create location descriptor fixture for BRAF catalytic domain"""
    params = {
        "id": "fusor.location_descriptor:NP_004324.2",
        "type": "LocationDescriptor",
        "location": {
            "sequence_id": "ga4gh:SQ.cQvw4UsHHRRlogxbWCB8W-mKD4AraM9y",
            "type": "SequenceLocation",
            "interval": {
                "start": {"type": "Number", "value": 458},
                "end": {
                    "type": "Number",
                    "value": 712,
                },
            },
        },
    }
    return LocationDescriptor(**params)


@pytest.fixture(scope="module")
def functional_domain_min(braf_gene_descr_min, location_descriptor_braf_domain):
    """Create functional domain test fixture."""
    params = {
        "status": "preserved",
        "label": "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "id": "interpro:IPR001245",
        "associated_gene": braf_gene_descr_min,
        "sequence_location": location_descriptor_braf_domain,
    }
    return FunctionalDomain(**params)


@pytest.fixture(scope="module")
def functional_domain(braf_gene_descr, location_descriptor_braf_domain):
    """Create functional domain test fixture."""
    params = {
        "status": "preserved",
        "label": "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "id": "interpro:IPR001245",
        "associated_gene": braf_gene_descr,
        "sequence_location": location_descriptor_braf_domain,
    }
    return FunctionalDomain(**params)


@pytest.fixture(scope="module")
def functional_domain_seq_id(
    braf_gene_descr_min, location_descriptor_braf_domain_seq_id
):
    """Create functional domain test fixture."""
    params = {
        "status": "preserved",
        "label": "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "id": "interpro:IPR001245",
        "associated_gene": braf_gene_descr_min,
        "sequence_location": location_descriptor_braf_domain_seq_id,
    }
    return FunctionalDomain(**params)


@pytest.fixture(scope="module")
def regulatory_element(braf_gene_descr):
    """Create regulatory element test fixture."""
    params = {
        "type": "RegulatoryElement",
        "regulatory_class": "promoter",
        "associated_gene": braf_gene_descr,
    }
    return RegulatoryElement(**params)


@pytest.fixture(scope="module")
def regulatory_element_min(braf_gene_descr_min):
    """Create regulatory element test fixture with minimal gene descriptor."""
    params = {"regulatory_class": "promoter", "associated_gene": braf_gene_descr_min}
    return RegulatoryElement(**params)


@pytest.fixture(scope="module")
def location_descriptor_tpm3():
    """Create location descriptor test fixture."""
    params = {
        "id": "fusor.location_descriptor:NM_152263.3",
        "type": "LocationDescriptor",
        "location": {
            "sequence_id": "refseq:NM_152263.3",
            "type": "SequenceLocation",
            "interval": {
                "start": {"type": "Number", "value": 154170398},
                "end": {"type": "Number", "value": 154170399},
                "type": "SequenceInterval",
            },
        },
    }
    return LocationDescriptor(**params)


@pytest.fixture(scope="module")
def templated_sequence_element():
    """Create test fixture for templated sequence element"""
    params = {
        "type": "TemplatedSequenceElement",
        "region": {
            "id": "fusor.location_descriptor:NC_000001.11",
            "type": "LocationDescriptor",
            "location": {
                "type": "SequenceLocation",
                "sequence_id": "refseq:NC_000001.11",
                "interval": {
                    "type": "SequenceInterval",
                    "start": {"type": "Number", "value": 99},
                    "end": {"type": "Number", "value": 150},
                },
            },
        },
        "strand": "+",
    }
    return TemplatedSequenceElement(**params)


@pytest.fixture()
def templated_sequence_element_ensg():
    """Create test fixture using non-seqrepo-recognized sequence ID"""
    params = {
        "type": "TemplatedSequenceElement",
        "region": {
            "id": "fusor.location_descriptor:ENSG00000157764",
            "type": "LocationDescriptor",
            "location": {
                "type": "SequenceLocation",
                "sequence_id": "ensembl:ENSG00000157764",
                "interval": {
                    "type": "SequenceInterval",
                    "start": {"type": "Number", "value": 140719328},
                    "end": {"type": "Number", "value": 140719400},
                },
            },
        },
        "strand": "-",
    }
    return TemplatedSequenceElement(**params)


@pytest.fixture(scope="module")
def templated_sequence_element_custom_id():
    """Create test fixture using custom (ie unable to coerce namespace)
    sequence identifier.
    """
    params = {
        "type": "TemplatedSequenceElement",
        "region": {
            "id": "fusor.location_descriptor:custom_ID__1",
            "type": "LocationDescriptor",
            "location": {
                "type": "SequenceLocation",
                "sequence_id": "sequence.id:custom_ID__1",
                "interval": {
                    "type": "SequenceInterval",
                    "start": {"type": "Number", "value": 200},
                    "end": {"type": "Number", "value": 300},
                },
            },
        },
        "strand": "+",
    }
    return TemplatedSequenceElement(**params)


@pytest.fixture(scope="module")
def transcript_segment_element():
    """Create transcript segment element test fixture"""
    params = {
        "type": "TranscriptSegmentElement",
        "exon_end": 8,
        "exon_end_offset": 0,
        "exon_start": 1,
        "exon_start_offset": 0,
        "gene_descriptor": {
            "gene_id": "hgnc:12012",
            "id": "normalize.gene:TPM3",
            "label": "TPM3",
            "type": "GeneDescriptor",
        },
        "transcript": "refseq:NM_152263.3",
        "element_genomic_end": {
            "id": "fusor.location_descriptor:NC_000001.11",
            "label": "NC_000001.11",
            "location": {
                "interval": {
                    "end": {"type": "Number", "value": 154170400},
                    "start": {"type": "Number", "value": 154170399},
                    "type": "SequenceInterval",
                },
                "sequence_id": "refseq:NC_000001.11",
                "type": "SequenceLocation",
            },
            "type": "LocationDescriptor",
        },
        "element_genomic_start": {
            "id": "fusor.location_descriptor:NC_000001.11",
            "label": "NC_000001.11",
            "location": {
                "interval": {
                    "end": {"type": "Number", "value": 154192136},
                    "start": {"type": "Number", "value": 154192135},
                    "type": "SequenceInterval",
                },
                "sequence_id": "refseq:NC_000001.11",
                "type": "SequenceLocation",
            },
            "type": "LocationDescriptor",
        },
    }
    return TranscriptSegmentElement(**params)


@pytest.fixture(scope="module")
def mane_transcript_segment_element():
    """Create transcript segment element test fixture"""
    params = {
        "type": "TranscriptSegmentElement",
        "exon_end": None,
        "exon_end_offset": None,
        "exon_start": 2,
        "exon_start_offset": 0,
        "gene_descriptor": {
            "gene_id": "hgnc:12761",
            "id": "normalize.gene:WEE1",
            "label": "WEE1",
            "type": "GeneDescriptor",
        },
        "transcript": "refseq:NM_003390.4",
        "element_genomic_end": None,
        "element_genomic_start": {
            "id": "fusor.location_descriptor:NC_000011.10",
            "label": "NC_000011.10",
            "location": {
                "interval": {
                    "end": {"type": "Number", "value": 9576094},
                    "start": {"type": "Number", "value": 9576093},
                    "type": "SequenceInterval",
                },
                "sequence_id": "refseq:NC_000011.10",
                "type": "SequenceLocation",
            },
            "type": "LocationDescriptor",
        },
    }
    return TranscriptSegmentElement(**params)


@pytest.fixture()
def fusion_ensg_sequence_id(templated_sequence_element_ensg):
    """Create fixture using Ensemble gene ID."""
    params = {
        "type": "CategoricalFusion",
        "structural_elements": [
            templated_sequence_element_ensg,
            {"type": "MultiplePossibleGenesElement"},
        ],
        "r_frame_preserved": True,
        "regulatory_element": None,
    }
    return CategoricalFusion(**params)


def compare_gene_descriptor(actual: dict, expected: dict):
    """Test that actual and expected gene descriptors match."""
    assert actual["id"] == expected["id"]
    assert actual["type"] == expected["type"]
    assert actual["gene_id"] == expected["gene_id"]
    assert actual["label"] == expected["label"]
    if expected["xrefs"]:
        assert set(actual["xrefs"]) == set(expected["xrefs"]), "xrefs"
    else:
        assert actual["xrefs"] == expected["xrefs"]
    if expected["alternate_labels"]:
        assert set(actual["alternate_labels"]) == set(
            expected["alternate_labels"]
        ), "alt labels"
    else:
        assert actual["alternate_labels"] == expected["alternate_labels"]
    assert "extensions" in actual
    if expected["extensions"]:
        assert len(actual["extensions"]) == len(
            expected["extensions"]
        ), "len of extensions"
        n_ext_correct = 0
        for expected_ext in expected["extensions"]:
            for actual_ext in actual["extensions"]:
                if actual_ext["name"] == expected_ext["name"]:
                    assert isinstance(actual_ext["value"], type(expected_ext["value"]))
                    if isinstance(expected_ext["value"], list) and not isinstance(
                        expected_ext["value"][0], dict
                    ):
                        assert set(actual_ext["value"]) == set(
                            expected_ext["value"]
                        ), f"{expected_ext['value']} value"
                    else:
                        assert actual_ext["value"] == expected_ext["value"]
                    assert actual_ext["type"] == expected_ext["type"]
                    n_ext_correct += 1
        assert n_ext_correct == len(
            expected["extensions"]
        ), "number of correct extensions"


def test_add_additional_fields(fusor_instance, fusion_example, fusion_ensg_sequence_id):
    """Test that add_additional_fields method works correctly."""
    fusion = CategoricalFusion(**fusion_example)

    expected_fusion = copy.deepcopy(fusion)
    expected_fusion.critical_functional_domains[
        0
    ].sequence_location.location_id = "ga4gh:VSL.2CWYzSpOJfZq7KW4VIUKeP5SJtepRar0"
    expected_fusion.critical_functional_domains[
        0
    ].sequence_location.location.sequence_id = (
        "ga4gh:SQ.q9CnK-HKWh9eqhOi8FlzR7M0pCmUrWPs"
    )
    expected_fusion.structural_elements[
        0
    ].element_genomic_start.location_id = "ga4gh:VSL.H0IOyJ-DB4jTbbSBjQFvuPvMrZHAWSrW"
    expected_fusion.structural_elements[
        0
    ].element_genomic_start.location.sequence_id = (
        "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"
    )
    expected_fusion.structural_elements[
        0
    ].element_genomic_end.location_id = "ga4gh:VSL.aarSLdMOQ8LoooPB2EoSth41yG_qRmDq"
    expected_fusion.structural_elements[
        0
    ].element_genomic_end.location.sequence_id = (
        "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"
    )
    expected_fusion.structural_elements[
        3
    ].region.location_id = "ga4gh:VSL.zd12pX_ju2gLq9a9UOYgM8AtbkuhnyUu"
    expected_fusion.structural_elements[
        3
    ].region.location.sequence_id = "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP"

    actual_fusion = fusor_instance.add_additional_fields(fusion)
    assert actual_fusion.model_dump() == expected_fusion.model_dump()

    # test handling of unrecognized sequence IDs
    expected_fusion = copy.deepcopy(fusion_ensg_sequence_id)
    fusion = fusor_instance.add_additional_fields(fusion_ensg_sequence_id)
    ts_reg = fusion.structural_elements[0].region
    assert ts_reg.location.sequence_id == "ensembl:ENSG00000157764"
    assert ts_reg.location_id == "ga4gh:VSL.dUll0TA05efQf0TsmcP03mtdGcpP9jPH"


def test_add_translated_sequence_id(fusor_instance, fusion_example):
    """Test that add_translated_sequence_id method works correctly."""
    fusion = CategoricalFusion(**fusion_example)

    expected_fusion = copy.deepcopy(fusion)
    expected_fusion.critical_functional_domains[
        0
    ].sequence_location.location.sequence_id = (
        "ga4gh:SQ.q9CnK-HKWh9eqhOi8FlzR7M0pCmUrWPs"
    )
    expected_fusion.structural_elements[
        0
    ].element_genomic_start.location.sequence_id = (
        "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"
    )
    expected_fusion.structural_elements[
        0
    ].element_genomic_end.location.sequence_id = (
        "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"
    )
    expected_fusion.structural_elements[
        3
    ].region.location.sequence_id = "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP"

    actual_fusion = fusor_instance.add_translated_sequence_id(fusion)
    assert actual_fusion.model_dump() == expected_fusion.model_dump()


def test_add_location_id(fusor_instance, fusion_example, exhaustive_example):
    """Test that add_location_id method works correctly."""
    fusion = fusor_instance.add_location_id(CategoricalFusion(**fusion_example))
    actual = CategoricalFusion(**exhaustive_example)

    assert (
        fusion.critical_functional_domains[0].sequence_location.location_id
        == actual.critical_functional_domains[0].sequence_location.location_id
    )
    assert (
        fusion.structural_elements[0].element_genomic_start.location_id
        == actual.structural_elements[0].element_genomic_start.location_id
    )
    assert (
        fusion.structural_elements[0].element_genomic_end.location_id
        == actual.structural_elements[0].element_genomic_end.location_id
    )
    assert (
        fusion.structural_elements[3].region.location_id
        == actual.structural_elements[3].region.location_id
    )


def test__normalized_gene_descriptor(fusor_instance):
    """Test that _normalized_gene_descriptor works correctly."""
    # Actual response is tested in test_add_gene_descriptor
    resp = fusor_instance._normalized_gene_descriptor("BRAF")
    assert resp[0]
    assert resp[1] is None
    assert isinstance(resp[0], GeneDescriptor)

    resp = fusor_instance._normalized_gene_descriptor("B R A F")
    assert resp[0] is None
    assert resp[1] == "gene-normalizer unable to normalize B R A F"


def test_add_gene_descriptor(fusor_instance, exhaustive_example, fusion_example):
    """Test that add_gene_descriptor method works correctly."""
    expected_fusion = CategoricalFusion(**exhaustive_example)
    actual = CategoricalFusion(**fusion_example)
    fusor_instance.add_translated_sequence_id(actual)
    fusor_instance.add_location_id(actual)
    fusor_instance.add_gene_descriptor(actual)

    e_gds = set()
    t_gds = set()
    for e_field in [
        expected_fusion.critical_functional_domains,
        expected_fusion.structural_elements,
    ]:
        for t_field in [actual.critical_functional_domains, actual.structural_elements]:
            for e_obj in e_field:
                for t_obj in t_field:
                    if "gene_descriptor" in e_obj.model_fields:
                        e_gd = e_obj.gene_descriptor.label
                        e_gds.add(e_gd)
                        if "gene_descriptor" in t_obj.model_fields:
                            t_gd = t_obj.gene_descriptor.label
                            t_gds.add(t_gd)
                            if e_gd == t_gd:
                                compare_gene_descriptor(
                                    t_obj.gene_descriptor.model_dump(),
                                    e_obj.gene_descriptor.model_dump(),
                                )
    assert t_gds == e_gds

    compare_gene_descriptor(
        actual.regulatory_element.associated_gene.model_dump(),
        expected_fusion.regulatory_element.associated_gene.model_dump(),
    )


def test_fusion(
    fusor_instance,
    linker_element,
    templated_sequence_element,
    transcript_segment_element,
    functional_domain,
):
    """Test that fusion methods work correctly."""
    # infer type from properties
    f = fusor_instance.fusion(
        structural_elements=[
            templated_sequence_element,
            linker_element,
            UnknownGeneElement(),
        ],
        causative_event={
            "type": "CausativeEvent",
            "event_type": "rearrangement",
            "event_description": "chr2:g.pter_8,247,756::chr11:g.15,825,273_cen_qter (der11) and chr11:g.pter_15,825,272::chr2:g.8,247,757_cen_qter (der2)",
        },
        assay={
            "type": "Assay",
            "method_uri": "pmid:33576979",
            "assay_id": "obi:OBI_0003094",
            "assay_name": "fluorescence in-situ hybridization assay",
            "fusion_detection": "inferred",
        },
    )
    assert isinstance(f, AssayedFusion)
    f = fusor_instance.fusion(
        structural_elements=[
            transcript_segment_element,
            MultiplePossibleGenesElement(),
        ],
        critical_functional_domains=[functional_domain],
    )
    assert isinstance(f, CategoricalFusion)

    # catch conflicting property args
    with pytest.raises(FUSORParametersException) as excinfo:
        f = fusor_instance.fusion(
            structural_elements=[
                transcript_segment_element,
                UnknownGeneElement(),
            ],
            causative_event="rearrangement",
            critical_functional_domains=[functional_domain],
        )
    assert str(excinfo.value) == "Received conflicting attributes"

    # handle indeterminate type
    with pytest.raises(FUSORParametersException) as excinfo:
        f = fusor_instance.fusion(
            structural_elements=[
                transcript_segment_element,
                templated_sequence_element,
            ]
        )
    assert str(excinfo.value) == "Unable to determine fusion type"

    # handle both type parameter options
    f = fusor_instance.fusion(
        fusion_type="AssayedFusion",
        structural_elements=[
            templated_sequence_element,
            linker_element,
            UnknownGeneElement(),
        ],
        causative_event={
            "type": "CausativeEvent",
            "event_type": "rearrangement",
        },
        assay={
            "type": "Assay",
            "method_uri": "pmid:33576979",
            "assay_id": "obi:OBI_0003094",
            "assay_name": "fluorescence in-situ hybridization assay",
            "fusion_detection": "inferred",
        },
    )
    assert isinstance(f, AssayedFusion)
    f = fusor_instance.fusion(
        type="CategoricalFusion",
        structural_elements=[
            transcript_segment_element,
            MultiplePossibleGenesElement(),
        ],
        critical_functional_domains=[functional_domain],
    )
    assert isinstance(f, CategoricalFusion)

    # catch and pass on validation errors
    with pytest.raises(FUSORParametersException) as excinfo:
        f = fusor_instance.fusion(
            fusion_type="CategoricalFusion", structural_elements=[linker_element]
        )
    msg = "Fusions must contain >= 2 structural elements, or >=1 structural element and a regulatory element"
    assert msg in str(excinfo.value)


@pytest.mark.asyncio()
async def test_transcript_segment_element(
    fusor_instance, transcript_segment_element, mane_transcript_segment_element
):
    """Test that transcript_segment_element method works correctly"""
    # Transcript Input
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3", exon_start=1, exon_end=8, tx_to_genomic_coords=True
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == transcript_segment_element.model_dump()

    # Genomic input, residue
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3",
        start=154192136,
        end=154170399,
        chromosome="NC_000001.11",
        tx_to_genomic_coords=False,
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == transcript_segment_element.model_dump()

    # Genomic input, inter-residue
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3",
        start=154192135,
        end=154170399,
        chromosome="NC_000001.11",
        tx_to_genomic_coords=False,
        residue_mode="inter-residue",
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == transcript_segment_element.model_dump()

    # Transcript Input
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3",
        exon_start=1,
        exon_end=8,
        gene="TPM3",
        tx_to_genomic_coords=True,
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == transcript_segment_element.model_dump()

    expected = copy.deepcopy(transcript_segment_element)
    expected.element_genomic_start.location.sequence_id = (
        "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"
    )
    expected.element_genomic_end.location.sequence_id = (
        expected.element_genomic_start.location.sequence_id
    )

    # Transcript Input
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3",
        exon_start=1,
        exon_end=8,
        tx_to_genomic_coords=True,
        seq_id_target_namespace="ga4gh",
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == expected.model_dump()

    # Genomic input
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3",
        start=154192136,
        end=154170399,
        chromosome="NC_000001.11",
        tx_to_genomic_coords=False,
        seq_id_target_namespace="ga4gh",
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == expected.model_dump()

    expected.exon_end_offset = -5
    expected.element_genomic_end.location.interval.start.value = 154170404
    expected.element_genomic_end.location.interval.end.value = 154170405

    # Transcript Input
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3",
        exon_start=1,
        exon_end=8,
        exon_end_offset=-5,
        tx_to_genomic_coords=True,
        seq_id_target_namespace="ga4gh",
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == expected.model_dump()

    # Genomic Input
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3",
        start=154192136,
        end=154170404,
        chromosome="NC_000001.11",
        tx_to_genomic_coords=False,
        seq_id_target_namespace="ga4gh",
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == expected.model_dump()

    expected.exon_end = None
    expected.exon_end_offset = None
    expected.element_genomic_end = None

    # Transcript Input
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3",
        exon_start=1,
        tx_to_genomic_coords=True,
        seq_id_target_namespace="ga4gh",
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == expected.model_dump()

    # Genomic Input
    tsg = await fusor_instance.transcript_segment_element(
        transcript="NM_152263.3",
        start=154192136,
        chromosome="NC_000001.11",
        tx_to_genomic_coords=False,
        seq_id_target_namespace="ga4gh",
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == expected.model_dump()

    # MANE
    tsg = await fusor_instance.transcript_segment_element(
        tx_to_genomic_coords=False,
        chromosome="NC_000011.10",
        start=9576094,
        gene="WEE1",
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].model_dump() == mane_transcript_segment_element.model_dump()


def test_gene_element(fusor_instance, braf_gene_descr_min, braf_gene_descr):
    """Test that gene_element works correctly."""
    gc = fusor_instance.gene_element("BRAF", use_minimal_gene_descr=True)
    assert gc[0]
    assert gc[1] is None
    assert isinstance(gc[0], GeneElement)
    compare_gene_descriptor(
        gc[0].gene_descriptor.model_dump(), braf_gene_descr_min.model_dump()
    )

    gc = fusor_instance.gene_element("BRAF", use_minimal_gene_descr=False)
    assert gc[0]
    assert gc[1] is None
    assert isinstance(gc[0], GeneElement)
    compare_gene_descriptor(
        gc[0].gene_descriptor.model_dump(), braf_gene_descr.model_dump()
    )

    gc = fusor_instance.gene_element("BRA F", use_minimal_gene_descr=True)
    assert gc[0] is None
    assert gc[1] == "gene-normalizer unable to normalize BRA F"


def test_templated_sequence_element(
    fusor_instance,
    templated_sequence_element,
    templated_sequence_element_ensg,
    templated_sequence_element_custom_id,
):
    """Test that templated sequence element works correctly"""
    tsg = fusor_instance.templated_sequence_element(
        100, 150, "NC_000001.11", "+", residue_mode="residue"
    )
    assert tsg.model_dump() == templated_sequence_element.model_dump()

    tsg = fusor_instance.templated_sequence_element(
        99, 150, "NC_000001.11", "+", residue_mode="inter-residue"
    )
    assert tsg.model_dump() == templated_sequence_element.model_dump()

    expected = copy.deepcopy(templated_sequence_element.model_dump())
    expected["region"]["location"][
        "sequence_id"
    ] = "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"
    expected["region"]["location_id"] = "ga4gh:VSL.bL1N-PQfp4dGlEz6PEd34fGxdxo82Zkb"
    tsg = fusor_instance.templated_sequence_element(
        100,
        150,
        "NC_000001.11",
        "+",
        add_location_id=True,
        seq_id_target_namespace="ga4gh",
    )
    assert tsg.model_dump() == expected

    tsg = fusor_instance.templated_sequence_element(
        140719329, 140719400, "ENSG00000157764", "-"
    )
    assert tsg.model_dump() == templated_sequence_element_ensg.model_dump()

    # test untranslateable sequence ID
    # adds "ensembl" namespace but unable to translate to ga4gh digest ID
    expected = copy.deepcopy(templated_sequence_element_ensg.model_dump())
    tsg = fusor_instance.templated_sequence_element(
        140719329, 140719400, "ENSG00000157764", "-", seq_id_target_namespace="ga4gh"
    )
    assert tsg.model_dump() == expected

    # test in-house/bespoke sequence ID
    # can't coerce namespace or translate to ga4gh ID
    expected = copy.deepcopy(templated_sequence_element_custom_id.model_dump())
    tsg = fusor_instance.templated_sequence_element(
        200,
        300,
        "custom_ID__1",
        "+",
        residue_mode="inter-residue",
        seq_id_target_namespace="ga4gh",
    )
    assert tsg.model_dump() == expected


def test_linker_element(fusor_instance, linker_element):
    """Test that linker_element method works correctly."""
    lc = fusor_instance.linker_element("act")
    assert lc[0]
    assert lc[1] is None
    assert lc[0].model_dump() == linker_element.model_dump()

    lc = fusor_instance.linker_element("bob!")
    assert lc[0] is None
    assert "String should match pattern '^[A-Z*\\-]*$'" in lc[1]


def test_unknown_gene_element(fusor_instance):
    """Test that unknown_gene_element method works correctly."""
    unknown_gc = fusor_instance.unknown_gene_element()
    assert unknown_gc.model_dump() == UnknownGeneElement().model_dump()


def test_multiple_possible_genes_element(fusor_instance):
    """Test that test_multiple_possible_genes_element method works correctly."""
    mult_gene = fusor_instance.multiple_possible_genes_element()
    assert mult_gene.model_dump() == MultiplePossibleGenesElement().model_dump()


def test_functional_domain(
    fusor_instance, functional_domain, functional_domain_min, functional_domain_seq_id
):
    """Test that functional_domain method works correctly"""

    def compare_domains(actual, expected):
        """Compare actual and expected functional domain data"""
        assert actual[0]
        assert actual[1] is None
        actual = actual[0].model_dump()
        expected = expected.model_dump()
        assert actual.keys() == expected.keys()
        for key in expected:
            if key == "associated_gene":
                compare_gene_descriptor(actual[key], expected[key])
            elif key == "sequence_location":
                act_ld = actual["sequence_location"]
                exp_ld = expected["sequence_location"]
                assert act_ld["id"] == exp_ld["id"]
                assert act_ld["type"] == exp_ld["type"]
                assert act_ld["location"]["type"] == exp_ld["location"]["type"]
                assert (
                    act_ld["location"]["sequence_id"]
                    == exp_ld["location"]["sequence_id"]
                )
                act_int = act_ld["location"]["interval"]
                exp_int = exp_ld["location"]["interval"]
                assert act_int["type"] == exp_int["type"]
                assert act_int["start"]["type"] == exp_int["start"]["type"]
                assert act_int["start"]["value"] == exp_int["start"]["value"]
                assert act_int["end"]["type"] == exp_int["end"]["type"]
                assert act_int["end"]["value"] == exp_int["end"]["value"]
            else:
                assert actual[key] == expected[key]

    cd = fusor_instance.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245",
        "BRAF",
        "NP_004324.2",
        458,
        712,
        use_minimal_gene_descr=False,
    )
    compare_domains(cd, functional_domain)

    cd = fusor_instance.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245",
        "BRAF",
        "NP_004324.2",
        458,
        712,
        use_minimal_gene_descr=True,
    )
    compare_domains(cd, functional_domain_min)

    cd = fusor_instance.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245",
        "BRAF",
        "NP_004324.2",
        458,
        712,
        seq_id_target_namespace="ga4gh",
        use_minimal_gene_descr=True,
    )
    compare_domains(cd, functional_domain_seq_id)

    cd = fusor_instance.functional_domain(
        "preserveded",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245",
        "BRAF",
        "NP_004324.2",
        458,
        712,
        seq_id_target_namespace="ga4gh",
        use_minimal_gene_descr=True,
    )
    assert cd[0] is None
    assert "Input should be 'lost' or 'preserved'" in cd[1]

    # check for protein accession
    cd = fusor_instance.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245",
        "BRAF",
        "NM_004333.4",
        458,
        712,
        seq_id_target_namespace="ga4gh",
        use_minimal_gene_descr=True,
    )
    assert cd[0] is None
    assert "Sequence_id must be a protein accession." in cd[1]

    # check for recognized protein accession
    accession = "NP_9999.999"
    cd = fusor_instance.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245",
        "BRAF",
        accession,
        458,
        712,
        seq_id_target_namespace="ga4gh",
        use_minimal_gene_descr=True,
    )
    assert cd[0] is None
    assert f"Accession, {accession}, not found in SeqRepo" in cd[1]

    # check that coordinates exist on sequence
    cd = fusor_instance.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245",
        "BRAF",
        "NP_004324.2",
        458,
        712000,
        seq_id_target_namespace="ga4gh",
        use_minimal_gene_descr=True,
    )
    assert cd[0] is None
    assert (
        "End inter-residue coordinate (711999) is out of index on "
        "NP_004324.2" in cd[1]
    )


def test_regulatory_element(fusor_instance, regulatory_element, regulatory_element_min):
    """Test regulatory_element method."""

    def compare_re(actual, expected):
        """Compare actual and expected regulatory element results."""
        assert actual[0]
        assert actual[1] is None
        actual = actual[0].model_dump()
        expected = expected.model_dump()
        assert actual.keys() == expected.keys()
        assert actual["type"] == expected["type"]
        compare_gene_descriptor(actual["associated_gene"], expected["associated_gene"])

    re = fusor_instance.regulatory_element(RegulatoryClass.PROMOTER, "BRAF")
    compare_re(re, regulatory_element_min)

    re = fusor_instance.regulatory_element(RegulatoryClass.PROMOTER, "BRAF", False)
    compare_re(re, regulatory_element)


def test__location_descriptor(fusor_instance, location_descriptor_tpm3):
    """Test that _location_descriptor method works correctly."""
    ld = fusor_instance._location_descriptor(154170398, 154170399, "NM_152263.3")
    assert ld.model_dump() == location_descriptor_tpm3.model_dump()

    expected = copy.deepcopy(location_descriptor_tpm3)
    expected.location.sequence_id = "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT"
    ld = fusor_instance._location_descriptor(
        154170398, 154170399, "NM_152263.3", seq_id_target_namespace="ga4gh"
    )
    assert ld.model_dump() == expected.model_dump()

    expected.id = "ga4gh:VSL._1bRdL4I6EtpBvVK5RUaXb0NN3k0gpqa"
    ld = fusor_instance._location_descriptor(
        154170398,
        154170399,
        "NM_152263.3",
        seq_id_target_namespace="ga4gh",
        use_location_id=True,
    )
    assert ld.model_dump() == expected.model_dump()

    expected.location.sequence_id = "refseq:NM_152263.3"
    expected.id = "fusor.location_descriptor:refseq%3ANM_152263.3"
    ld = fusor_instance._location_descriptor(154170398, 154170399, "refseq:NM_152263.3")
    assert ld.model_dump() == expected.model_dump()

    expected.id = "fusor.location_descriptor:example_label"
    expected.label = "example_label"
    ld = fusor_instance._location_descriptor(
        154170398, 154170399, "refseq:NM_152263.3", label="example_label"
    )
    assert ld.model_dump() == expected.model_dump()
