"""Module for testing the FUSOR class."""
import copy

import pytest
from ga4gh.vrsatile.pydantic.vrsatile_models import GeneDescriptor, \
    LocationDescriptor

from fusor import FUSOR
from fusor.exceptions import IDTranslationException
from fusor.models import Fusion, TemplatedSequenceComponent, \
    TranscriptSegmentComponent, LinkerComponent, UnknownGeneComponent, \
    AnyGeneComponent, FunctionalDomain, GeneComponent, RegulatoryElement, \
    RegulatoryElementType


@pytest.fixture(scope="module")
def fusor():
    """Create test fixture for fusor object"""
    return FUSOR()


@pytest.fixture(scope="module")
def braf_gene_descr_min():
    """Create minimal gene descriptor for BRAF"""
    return GeneDescriptor(
        id="normalize.gene:BRAF",
        label="BRAF",
        gene_id="hgnc:1097"
    )


@pytest.fixture(scope="module")
def braf_gene_descr():
    """Create gene descriptor for BRAF."""
    params = {
        "id": "normalize.gene:BRAF",
        "type": "GeneDescriptor",
        "gene_id": "hgnc:1097",
        "label": "BRAF",
        "xrefs": [
            "ncbigene:673",
            "ensembl:ENSG00000157764"
        ],
        "alternate_labels": [
            "BRAF1",
            "RAFB1",
            "B-raf",
            "NS7",
            "B-RAF1"
        ],
        "extensions": [
            {
                "name": "approved_name",
                "value": "B-Raf proto-oncogene, serine/threonine kinase",
                "type": "Extension"
            },
            {
                "name": "symbol_status",
                "value": "approved",
                "type": "Extension"
            },
            {
                "name": "associated_with",
                "value": [
                    "ccds:CCDS5863",
                    "iuphar:1943",
                    "orphanet:119066",
                    "cosmic:BRAF",
                    "pubmed:2284096",
                    "ucsc:uc003vwc.5",
                    "omim:164757",
                    "refseq:NM_004333",
                    "ccds:CCDS87555",
                    "uniprot:P15056",
                    "ena.embl:M95712",
                    "vega:OTTHUMG00000157457",
                    "pubmed:1565476"
                ],
                "type": "Extension"
            },
            {
                "name": "chromosome_location",
                "value": {
                    "_id": "ga4gh:VCL.O6yCQ1cnThOrTfK9YUgMlTfM6HTqbrKw",
                    "type": "ChromosomeLocation",
                    "species_id": "taxonomy:9606",
                    "chr": "7",
                    "interval": {
                        "end": "q34",
                        "start": "q34",
                        "type": "CytobandInterval"
                    }
                },
                "type": "Extension"
            }
        ]
    }
    return GeneDescriptor(**params)


@pytest.fixture(scope="module")
def linker_component():
    """Create linker component test fixture."""
    params = {
        "linker_sequence": {
            "id": "fusor.sequence:ACT",
            "sequence": "ACT",
            "residue_type": "SO:0000348",
            "type": "SequenceDescriptor"
        },
        "component_type": "linker_sequence"
    }
    return LinkerComponent(**params)


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
                "start": {
                    "type": "Number",
                    "value": 458
                },
                "end": {
                    "type": "Number",
                    "value": 712,
                }
            }
        }
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
                "start": {
                    "type": "Number",
                    "value": 458
                },
                "end": {
                    "type": "Number",
                    "value": 712,
                }
            }
        }
    }
    return LocationDescriptor(**params)


@pytest.fixture(scope="module")
def functional_domain_min(braf_gene_descr_min,
                          location_descriptor_braf_domain):
    """Create functional domain test fixture."""
    params = {
        "status": "preserved",
        "name": "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "id": "interpro:IPR001245",
        "gene_descriptor": braf_gene_descr_min,
        "location_descriptor": location_descriptor_braf_domain
    }
    return FunctionalDomain(**params)


@pytest.fixture(scope="module")
def functional_domain(braf_gene_descr, location_descriptor_braf_domain):
    """Create functional domain test fixture."""
    params = {
        "status": "preserved",
        "name": "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "id": "interpro:IPR001245",
        "gene_descriptor": braf_gene_descr,
        "location_descriptor": location_descriptor_braf_domain
    }
    return FunctionalDomain(**params)


@pytest.fixture(scope="module")
def functional_domain_seq_id(braf_gene_descr_min,
                             location_descriptor_braf_domain_seq_id):
    """Create functional domain test fixture."""
    params = {
        "status": "preserved",
        "name": "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "id": "interpro:IPR001245",
        "gene_descriptor": braf_gene_descr_min,
        "location_descriptor": location_descriptor_braf_domain_seq_id
    }
    return FunctionalDomain(**params)


@pytest.fixture(scope="module")
def regulatory_element(braf_gene_descr):
    """Create regulatory element test fixture."""
    params = {
        "type": "promoter",
        "gene_descriptor": braf_gene_descr
    }
    return RegulatoryElement(**params)


@pytest.fixture(scope="module")
def regulatory_element_min(braf_gene_descr_min):
    """Create regulatory element test fixture with minimal gene descriptor."""
    params = {
        "type": "promoter",
        "gene_descriptor": braf_gene_descr_min
    }
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
    return LocationDescriptor(**params)


@pytest.fixture(scope="module")
def templated_sequence_component(location_descriptor_tpm3):
    """Create test fixture for templated sequence component"""
    params = {
        "component_type": "templated_sequence",
        "region": location_descriptor_tpm3.dict(exclude_none=True),
        "strand": "+"
    }
    return TemplatedSequenceComponent(**params)


@pytest.fixture(scope="module")
def transcript_segment_component():
    """Create transcript segment component test fixture"""
    params = {
        "component_type": "transcript_segment",
        "exon_end": 8,
        "exon_end_offset": 0,
        "exon_start": 1,
        "exon_start_offset": 0,
        "gene_descriptor": {
            "gene_id": "hgnc:12012",
            "id": "normalize.gene:TPM3",
            "label": "TPM3",
            "type": "GeneDescriptor"
        },
        "transcript": "refseq:NM_152263.3",
        "component_genomic_end": {
            "id": "fusor.location_descriptor:NC_000001.11",
            "label": "NC_000001.11",
            "location": {
                "interval": {
                    "end": {"type": "Number", "value": 154170400},
                    "start": {"type": "Number", "value": 154170399},
                    "type": "SequenceInterval"
                },
                "sequence_id": "refseq:NC_000001.11",
                "type": "SequenceLocation"
            },
            "type": "LocationDescriptor"
        },
        "component_genomic_start": {
            "id": "fusor.location_descriptor:NC_000001.11",
            "label": "NC_000001.11",
            "location": {
                "interval": {
                    "end": {"type": "Number", "value": 154192136},
                    "start": {"type": "Number", "value": 154192135},
                    "type": "SequenceInterval"
                },
                "sequence_id": "refseq:NC_000001.11",
                "type": "SequenceLocation"
            },
            "type": "LocationDescriptor"
        }
    }
    return TranscriptSegmentComponent(**params)


@pytest.fixture(scope="module")
def mane_transcript_segment_component():
    """Create transcript segment component test fixture"""
    params = {
        "component_type": "transcript_segment",
        "exon_end": None,
        "exon_end_offset": None,
        "exon_start": 2,
        "exon_start_offset": 0,
        "gene_descriptor": {
            "gene_id": "hgnc:12761",
            "id": "normalize.gene:WEE1",
            "label": "WEE1",
            "type": "GeneDescriptor"
        },
        "transcript": "refseq:NM_003390.4",
        "component_genomic_end": None,
        "component_genomic_start": {
            "id": "fusor.location_descriptor:NC_000011.10",
            "label": "NC_000011.10",
            "location": {
                "interval": {
                    "end": {"type": "Number", "value": 9576094},
                    "start": {"type": "Number", "value": 9576093},
                    "type": "SequenceInterval"
                },
                "sequence_id": "refseq:NC_000011.10",
                "type": "SequenceLocation"
            },
            "type": "LocationDescriptor"
        }
    }
    return TranscriptSegmentComponent(**params)


@pytest.fixture(scope="module")
def fusion():
    """Create fusion test fixture."""
    return {
        "r_frame_preserved": True,
        "functional_domains": [
            {
                "status": "lost",
                "name": "Tyrosine-protein kinase, catalytic domain",
                "id": "interpro:IPR020635",
                "gene_descriptor": {
                    "id": "gene:ALK",
                    "gene_id": "hgnc:1837",
                    "label": "ALK"
                },
                "location_descriptor": {
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
            },
            {
                "component_type": "any_gene"
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


@pytest.fixture(scope="module")
def fusion_custom_sequence_id():
    """Create fixture using custom/unrecognized sequence ID"""
    params = {
        "structural_components": [
            {
                "component_type": "templated_sequence",
                "strand": "+",
                "region": {
                    "id": "fusor.location_descriptor:NM_152263.3",
                    "type": "LocationDescriptor",
                    "location": {
                        "sequence_id": "refseq:chr12",
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
            {"component_type": "any_gene"}
        ],
        "r_frame_preserved": True,
        "functional_domains": [],
        "causative_event": None,
        "regulatory_elements": []
    }
    return Fusion(**params)


def compare_gene_descriptor(actual, expected):
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
        assert set(actual["alternate_labels"]) == set(expected["alternate_labels"]), "alt labels"  # noqa: E501
    else:
        assert actual["alternate_labels"] == expected["alternate_labels"]
    assert "extensions" in actual.keys() and "extensions" in expected.keys()
    if expected["extensions"]:
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


def test_add_additional_fields(fusor, fusion_example,
                               fusion_custom_sequence_id):
    """Test that add_additional_fields method works correctly."""
    fusion = Fusion(**fusion_example)

    expected_fusion = copy.deepcopy(fusion)
    expected_fusion.functional_domains[0].location_descriptor.location_id = "ga4gh:VSL.2CWYzSpOJfZq7KW4VIUKeP5SJtepRar0"  # type: ignore # noqa: E501
    expected_fusion.functional_domains[0].location_descriptor.location.sequence_id = "ga4gh:SQ.q9CnK-HKWh9eqhOi8FlzR7M0pCmUrWPs"  # type: ignore # noqa: E501
    expected_fusion.structural_components[0].component_genomic_start.location_id = "ga4gh:VSL.H0IOyJ-DB4jTbbSBjQFvuPvMrZHAWSrW"  # type: ignore # noqa: E501
    expected_fusion.structural_components[0].component_genomic_start.location.sequence_id = "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"  # type: ignore # noqa: E501
    expected_fusion.structural_components[0].component_genomic_end.location_id = "ga4gh:VSL.aarSLdMOQ8LoooPB2EoSth41yG_qRmDq"  # type: ignore # noqa: E501
    expected_fusion.structural_components[0].component_genomic_end.location.sequence_id = "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"  # type: ignore # noqa: E501
    expected_fusion.structural_components[3].region.location_id = "ga4gh:VSL.zd12pX_ju2gLq9a9UOYgM8AtbkuhnyUu"  # type: ignore # noqa: E501
    expected_fusion.structural_components[3].region.location.sequence_id = "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP"  # type: ignore # noqa: E501

    actual_fusion = fusor.add_additional_fields(fusion)
    assert actual_fusion.dict() == expected_fusion.dict()

    # test handling of custom/unrecognized sequence IDs
    expected_fusion = copy.deepcopy(fusion_custom_sequence_id)
    fusion = fusor.add_additional_fields(fusion_custom_sequence_id)
    ts_reg = fusion.structural_components[0].region
    assert ts_reg.location.sequence_id == "refseq:chr12"
    assert ts_reg.location_id == "ga4gh:VSL.K6sdndDkb0wSEZtTGbAqAJoEooe_2BA7"


def test_add_translated_sequence_id(fusor, fusion_example):
    """Test that add_translated_sequence_id method works correctly."""
    fusion = Fusion(**fusion_example)

    expected_fusion = copy.deepcopy(fusion)
    expected_fusion.functional_domains[0].location_descriptor.location.sequence_id = "ga4gh:SQ.q9CnK-HKWh9eqhOi8FlzR7M0pCmUrWPs"  # type: ignore # noqa: E501
    expected_fusion.structural_components[0].component_genomic_start.location.sequence_id = "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"  # type: ignore # noqa: E501
    expected_fusion.structural_components[0].component_genomic_end.location.sequence_id = "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"  # type: ignore # noqa: E501
    expected_fusion.structural_components[3].region.location.sequence_id = "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP"  # type: ignore # noqa: E501

    actual_fusion = fusor.add_translated_sequence_id(fusion)
    assert actual_fusion.dict() == expected_fusion.dict()


def test_add_location_id(fusor, fusion_example):
    """Test that add_location_id method works correctly."""
    fusion = Fusion(**fusion_example)

    expected_fusion = copy.deepcopy(fusion)
    expected_fusion.functional_domains[0].location_descriptor.location_id = "ga4gh:VSL.hQKhk6ZOOYZAmShXrzhfb6H3j65ovsKu"  # type: ignore # noqa: E501
    expected_fusion.structural_components[0].component_genomic_start.location_id = "ga4gh:VSL.n7i6VMRAuSgAjwVopxhWAJdlPJMfk7KR"  # type: ignore # noqa: E501
    expected_fusion.structural_components[0].component_genomic_end.location_id = "ga4gh:VSL.wQ4TpNbsTPq_A-eQTL44gbP3f4fnp0vx"  # type: ignore # noqa: E501
    expected_fusion.structural_components[3].region.location_id = "ga4gh:VSL.eHfgOlEjzNRiYZBzeCvg7Ru9N-YxuFT-"  # type: ignore # noqa: E501

    actual_fusion = fusor.add_location_id(fusion)
    assert actual_fusion.dict() == expected_fusion.dict()


def test_translate_identifier(fusor):
    """Test that translate_identifier method works correctly."""
    expected = "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT"
    identifier = fusor.translate_identifier("NM_152263.3")
    assert identifier == expected

    identifier = fusor.translate_identifier("refseq:NM_152263.3")
    assert identifier == expected

    # test non-default target
    identifier = fusor.translate_identifier(
        "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT",
        "refseq"
    )
    assert identifier == "refseq:NM_152263.3"

    # test no namespace
    with pytest.raises(IDTranslationException):
        identifier = fusor.translate_identifier("152263.3")

    # test unrecognized namespace
    with pytest.raises(IDTranslationException):
        identifier = fusor.translate_identifier("fake_namespace:NM_152263.3")


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
    actual = Fusion(**fusion)
    fusor.add_translated_sequence_id(actual)
    fusor.add_location_id(actual)
    fusor.add_gene_descriptor(actual)

    e_gds = set()
    t_gds = set()
    for e_field in [expected_fusion.functional_domains,
                    expected_fusion.structural_components,
                    expected_fusion.regulatory_elements]:
        for t_field in [actual.functional_domains,
                        actual.structural_components,
                        actual.regulatory_elements]:
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


def test_fusion(fusor, linker_component, templated_sequence_component,
                transcript_segment_component, functional_domain):
    """Test that fusion method works correctly."""
    f = fusor.fusion([templated_sequence_component,
                      linker_component, UnknownGeneComponent()])
    assert isinstance(f[0], Fusion)
    assert f[1] is None

    f = fusor.fusion([transcript_segment_component, AnyGeneComponent()],
                     functional_domains=[functional_domain])
    assert isinstance(f[0], Fusion)
    assert f[1] is None

    f = fusor.fusion([linker_component])
    assert f[0] is None
    assert "Fusion must contain at least 2 structural components" in f[1]


@pytest.mark.asyncio
async def test_transcript_segment_component(fusor,
                                            transcript_segment_component,
                                            mane_transcript_segment_component):
    """Test that transcript_segment_component method works correctly"""
    # Transcript Input
    tsg = await fusor.transcript_segment_component(
        transcript="NM_152263.3", exon_start=1, exon_end=8,
        tx_to_genomic_coords=True)
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].dict() == transcript_segment_component.dict()

    # Genomic input, residue
    tsg = await fusor.transcript_segment_component(
        transcript="NM_152263.3", start=154192136, end=154170399,
        chromosome="NC_000001.11", tx_to_genomic_coords=False
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].dict() == transcript_segment_component.dict()

    # Genomic input, inter-residue
    tsg = await fusor.transcript_segment_component(
        transcript="NM_152263.3", start=154192135, end=154170399,
        chromosome="NC_000001.11", tx_to_genomic_coords=False,
        residue_mode="inter-residue"
    )
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].dict() == transcript_segment_component.dict()

    # Transcript Input
    tsg = await fusor.transcript_segment_component(
        transcript="NM_152263.3", exon_start=1, exon_end=8, gene="TPM3",
        tx_to_genomic_coords=True)
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].dict() == transcript_segment_component.dict()

    expected = copy.deepcopy(transcript_segment_component)
    expected.component_genomic_start.location.sequence_id = \
        "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"
    expected.component_genomic_end.location.sequence_id = \
        expected.component_genomic_start.location.sequence_id

    # Transcript Input
    tsg = await fusor.transcript_segment_component(
        transcript="NM_152263.3", exon_start=1, exon_end=8,
        tx_to_genomic_coords=True, seq_id_target_namespace="ga4gh")
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].dict() == expected.dict()

    # Genomic input
    tsg = await fusor.transcript_segment_component(
        transcript="NM_152263.3", start=154192136, end=154170399,
        chromosome="NC_000001.11", tx_to_genomic_coords=False,
        seq_id_target_namespace="ga4gh")
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].dict() == expected.dict()

    expected.exon_end_offset = -5
    expected.component_genomic_end.location.interval.start.value = 154170404
    expected.component_genomic_end.location.interval.end.value = 154170405

    # Transcript Input
    tsg = await fusor.transcript_segment_component(
        transcript="NM_152263.3", exon_start=1, exon_end=8, exon_end_offset=-5,
        tx_to_genomic_coords=True, seq_id_target_namespace="ga4gh")
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].dict() == expected.dict()

    # Genomic Input
    tsg = await fusor.transcript_segment_component(
        transcript="NM_152263.3", start=154192136, end=154170404,
        chromosome="NC_000001.11", tx_to_genomic_coords=False,
        seq_id_target_namespace="ga4gh")
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].dict() == expected.dict()

    expected.exon_end = None
    expected.exon_end_offset = None
    expected.component_genomic_end = None

    # Transcript Input
    tsg = await fusor.transcript_segment_component(
        transcript="NM_152263.3", exon_start=1,
        tx_to_genomic_coords=True, seq_id_target_namespace="ga4gh")
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].dict() == expected.dict()

    # Genomic Input
    tsg = await fusor.transcript_segment_component(
        transcript="NM_152263.3", start=154192136,
        chromosome="NC_000001.11", tx_to_genomic_coords=False,
        seq_id_target_namespace="ga4gh")
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].dict() == expected.dict()

    # MANE
    tsg = await fusor.transcript_segment_component(
        tx_to_genomic_coords=False, chromosome="NC_000011.10",
        start=9576094, gene="WEE1")
    assert tsg[0]
    assert tsg[1] is None
    assert tsg[0].dict() == mane_transcript_segment_component.dict()


def test_gene_component(fusor, braf_gene_descr_min, braf_gene_descr):
    """Test that gene_component works correctly."""
    gc = fusor.gene_component("BRAF", use_minimal_gene_descr=True)
    assert gc[0]
    assert gc[1] is None
    assert isinstance(gc[0], GeneComponent)
    compare_gene_descriptor(gc[0].gene_descriptor.dict(),
                            braf_gene_descr_min.dict())

    gc = fusor.gene_component("BRAF", use_minimal_gene_descr=False)
    assert gc[0]
    assert gc[1] is None
    assert isinstance(gc[0], GeneComponent)
    compare_gene_descriptor(gc[0].gene_descriptor.dict(),
                            braf_gene_descr.dict())

    gc = fusor.gene_component("BRA F", use_minimal_gene_descr=True)
    assert gc[0] is None
    assert gc[1] == "gene-normalizer unable to normalize BRA F"


def test_templated_sequence_component(fusor, templated_sequence_component):
    """Test that templated sequence component works correctly"""
    tsg = fusor.templated_sequence_component(
        154170399, 154170399, "NM_152263.3", "+", residue_mode="residue"
    )
    assert tsg.dict() == templated_sequence_component.dict()

    tsg = fusor.templated_sequence_component(
        154170398, 154170399, "NM_152263.3", "+", residue_mode="inter-residue"
    )
    assert tsg.dict() == templated_sequence_component.dict()

    expected = copy.deepcopy(templated_sequence_component.dict())
    expected["region"]["location"]["sequence_id"] = "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT"  # noqa: E501
    expected["region"]["location_id"] = "ga4gh:VSL._1bRdL4I6EtpBvVK5RUaXb0NN3k0gpqa"  # noqa: E501
    tsg = fusor.templated_sequence_component(
        154170399, 154170399, "NM_152263.3", "+", add_location_id=True,
        seq_id_target_namespace="ga4gh"
    )
    assert tsg.dict() == expected


def test_linker_component(fusor, linker_component):
    """Test that linker_component method works correctly."""
    lc = fusor.linker_component("act")
    assert lc[0]
    assert lc[1] is None
    assert lc[0].dict() == linker_component.dict()

    lc = fusor.linker_component("bob!")
    assert lc[0] is None
    assert "sequence does not match regex '^[A-Za-z*\\-]*$'" in lc[1]


def test_unknown_gene_component(fusor):
    """Test that unknown_gene_component method works correctly."""
    unknown_gc = fusor.unknown_gene_component()
    assert unknown_gc.dict() == UnknownGeneComponent().dict()


def test_any_gene_component(fusor):
    """Test that any_gene component method works correctly."""
    any_gc = fusor.any_gene_component()
    assert any_gc.dict() == AnyGeneComponent().dict()


def test_functional_domain(fusor, functional_domain, functional_domain_min,
                           functional_domain_seq_id):
    """Test that functional_domain method works correctly"""

    def compare_domains(actual, expected):
        """Compare actual and expected functional domain data"""
        assert actual[0]
        assert actual[1] is None
        actual = actual[0].dict()
        expected = expected.dict()
        assert actual.keys() == expected.keys()
        for key in expected.keys():
            if key == "gene_descriptor":
                compare_gene_descriptor(actual[key], expected[key])
            elif key == "location_descriptor":
                act_ld = actual["location_descriptor"]
                exp_ld = expected["location_descriptor"]
                assert act_ld["id"] == exp_ld["id"]
                assert act_ld["type"] == exp_ld["type"]
                assert act_ld["location"]["type"] == exp_ld["location"]["type"]
                assert act_ld["location"]["sequence_id"] == \
                    exp_ld["location"]["sequence_id"]
                act_int = act_ld["location"]["interval"]
                exp_int = exp_ld["location"]["interval"]
                assert act_int["type"] == exp_int["type"]
                assert act_int["start"]["type"] == exp_int["start"]["type"]
                assert act_int["start"]["value"] == exp_int["start"]["value"]
                assert act_int["end"]["type"] == exp_int["end"]["type"]
                assert act_int["end"]["value"] == exp_int["end"]["value"]
            else:
                assert actual[key] == expected[key]

    cd = fusor.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245", "BRAF", "NP_004324.2", 458, 712,
        use_minimal_gene_descr=False)
    compare_domains(cd, functional_domain)

    cd = fusor.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245", "BRAF", "NP_004324.2", 458, 712,
        use_minimal_gene_descr=True)
    compare_domains(cd, functional_domain_min)

    cd = fusor.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245", "BRAF", "NP_004324.2", 458, 712,
        seq_id_target_namespace="ga4gh",
        use_minimal_gene_descr=True)
    compare_domains(cd, functional_domain_seq_id)

    cd = fusor.functional_domain(
        "preserveded",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245", "BRAF", "NP_004324.2", 458, 712,
        seq_id_target_namespace="ga4gh",
        use_minimal_gene_descr=True)
    assert cd[0] is None
    assert "value is not a valid enumeration member; permitted: " \
        "'lost', 'preserved'" in cd[1]

    # check for protein accession
    cd = fusor.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245", "BRAF", "NM_004333.4", 458, 712,
        seq_id_target_namespace="ga4gh",
        use_minimal_gene_descr=True)
    assert cd[0] is None
    assert "Sequence_id must be a protein accession." in cd[1]

    # check for recognized protein accession
    accession = "NP_9999.999"
    cd = fusor.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245", "BRAF",
        accession,
        458, 712,
        seq_id_target_namespace="ga4gh",
        use_minimal_gene_descr=True)
    assert cd[0] is None
    assert f"Accession, {accession}, not found in SeqRepo" in cd[1]

    # check that coordinates exist on sequence
    cd = fusor.functional_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245", "BRAF",
        "NP_004324.2",
        458, 712000,
        seq_id_target_namespace="ga4gh",
        use_minimal_gene_descr=True)
    assert cd[0] is None
    assert "End inter-residue coordinate (711999) is out of index on " \
        "NP_004324.2" in cd[1]


def test_regulatory_element(fusor, regulatory_element, regulatory_element_min):
    """Test regulatory_element method."""

    def compare_re(actual, expected):
        """Compare actual and expected regulatory element results."""
        assert actual[0]
        assert actual[1] is None
        actual = actual[0].dict()
        expected = expected.dict()
        assert actual.keys() == expected.keys()
        assert actual["type"] == expected["type"]
        compare_gene_descriptor(actual["gene_descriptor"],
                                expected["gene_descriptor"])

    re = fusor.regulatory_element(RegulatoryElementType.PROMOTER, "BRAF")
    compare_re(re, regulatory_element_min)

    re = fusor.regulatory_element(RegulatoryElementType.PROMOTER, "BRAF",
                                  False)
    compare_re(re, regulatory_element)


def test__location_descriptor(fusor, location_descriptor_tpm3):
    """Test that _location_descriptor method works correctly."""
    ld = fusor._location_descriptor(154170398, 154170399, "NM_152263.3")
    assert ld.dict() == location_descriptor_tpm3.dict()

    expected = copy.deepcopy(location_descriptor_tpm3)
    expected.location.sequence_id = "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT"
    ld = fusor._location_descriptor(154170398, 154170399, "NM_152263.3",
                                    seq_id_target_namespace="ga4gh")
    assert ld.dict() == expected.dict()

    expected.id = "ga4gh:VSL._1bRdL4I6EtpBvVK5RUaXb0NN3k0gpqa"
    ld = fusor._location_descriptor(154170398, 154170399, "NM_152263.3",
                                    seq_id_target_namespace="ga4gh",
                                    use_location_id=True)
    assert ld.dict() == expected.dict()

    expected.location.sequence_id = "refseq:NM_152263.3"
    expected.id = "fusor.location_descriptor:refseq%3ANM_152263.3"
    ld = fusor._location_descriptor(154170398, 154170399, "refseq:NM_152263.3")
    assert ld.dict() == expected.dict()

    expected.id = "fusor.location_descriptor:example_label"
    expected.label = "example_label"
    ld = fusor._location_descriptor(154170398, 154170399, "refseq:NM_152263.3",
                                    label="example_label")
    assert ld.dict() == expected.dict()


def test_generate_nomenclature(fusor, fusion):
    """Test that nomenclature generation is correct."""
    fusion_instance = Fusion(**fusion)
    nm = fusor.generate_nomenclature(fusion_instance)
    assert nm == "NM_152263.3(TPM3):e.1_8::ALK(hgnc:427)::ACGT::NC_000012.12:g.44908821_44908822(+)::?::*"  # noqa: E501
