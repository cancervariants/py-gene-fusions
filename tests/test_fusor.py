"""Module for testing the FUSOR class."""
import pytest
from ga4gh.vrsatile.pydantic.vrsatile_model import GeneDescriptor, \
    LocationDescriptor
from fusor import FUSOR
from fusor.models import Fusion, TemplatedSequenceComponent, \
    TranscriptSegmentComponent, LinkerComponent, UnknownGeneComponent, \
    CriticalDomain, GeneComponent
import copy


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
def critical_domain_min(braf_gene_descr_min):
    """Create critical domain test fixture."""
    params = {
        "status": "preserved",
        "name": "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "id": "interpro:IPR001245",
        "gene_descriptor": braf_gene_descr_min
    }
    return CriticalDomain(**params)


@pytest.fixture(scope="module")
def critical_domain(braf_gene_descr):
    """Create critical domain test fixture."""
    params = {
        "status": "preserved",
        "name": "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "id": "interpro:IPR001245",
        "gene_descriptor": braf_gene_descr
    }
    return CriticalDomain(**params)


@pytest.fixture(scope="module")
def location_descriptor():
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
def templated_sequence_component(location_descriptor):
    """Create test fixture for templated sequence component"""
    params = {
        "component_type": "templated_sequence",
        "region": location_descriptor.dict(exclude_none=True),
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
                    "end": {"type": "Number", "value": 154170399},
                    "start": {"type": "Number", "value": 154170398},
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


def test_add_additional_fields(fusor, fusion_example, fusion):
    """Test that add_additional_fields method works correctly."""
    expected_fusion = Fusion(**fusion_example)
    fusion = Fusion(**fusion)
    fusor.add_additional_fields(fusion)
    assert fusion.dict() == expected_fusion.dict()


def test_add_sequence_id(fusor, fusion_example, fusion):
    """Test that add_sequence_id method works correctly."""
    expected_fusion = Fusion(**fusion_example)
    actual = Fusion(**fusion)
    for structural_component in expected_fusion.structural_components:
        if isinstance(structural_component, TemplatedSequenceComponent):
            structural_component.region.location_id = None
        elif isinstance(structural_component, TranscriptSegmentComponent):
            if structural_component.component_genomic_start:
                structural_component.component_genomic_start.location_id = None
            if structural_component.component_genomic_end:
                structural_component.component_genomic_end.location_id = None
    fusor.add_sequence_id(actual)
    assert actual.dict() == expected_fusion.dict()


def test_add_location_id(fusor, fusion_example, fusion):
    """Test that add_location_id method works correctly."""
    expected_fusion = Fusion(**fusion_example)
    actual = Fusion(**fusion)
    fusor.add_sequence_id(actual)
    fusor.add_location_id(actual)
    assert actual.dict() == expected_fusion.dict()


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
    actual = Fusion(**fusion)
    fusor.add_sequence_id(actual)
    fusor.add_location_id(actual)
    fusor.add_gene_descriptor(actual)

    e_gds = set()
    t_gds = set()
    for e_field in [expected_fusion.protein_domains,
                    expected_fusion.structural_components,
                    expected_fusion.regulatory_elements]:
        for t_field in [actual.protein_domains,
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


def test_fusion(fusor, linker_component, templated_sequence_component):
    """Test that fusion method works correctly."""
    f = fusor.fusion([templated_sequence_component,
                      linker_component, UnknownGeneComponent()])
    assert isinstance(f[0], Fusion)
    assert f[1] is None

    f = fusor.fusion([linker_component])
    assert f[0] is None
    assert "Fusion must contain at least 2 structural components" in f[1]


@pytest.mark.asyncio
async def test_transcript_segment_component(fusor,
                                            transcript_segment_component):
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
    expected.component_genomic_end.location.interval.start.value = 154170403
    expected.component_genomic_end.location.interval.end.value = 154170404

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
    uknown_gc = fusor.unknown_gene_component()
    assert uknown_gc.dict() == UnknownGeneComponent().dict()


def test_critical_domain(fusor, critical_domain, critical_domain_min):
    """Test that critical_domain method works correctly"""

    def compare_cd(actual, expected):
        """Compare actual and expected critical domain data"""
        assert actual[0]
        assert actual[1] is None
        actual = actual[0].dict()
        expected = expected.dict()
        assert actual.keys() == expected.keys()
        for key in expected.keys():
            if key != "gene_descriptor":
                assert actual[key] == expected[key]
            else:
                compare_gene_descriptor(actual[key], expected[key])

    cd = fusor.critical_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245", "BRAF", use_minimal_gene_descr=True)
    compare_cd(cd, critical_domain_min)

    cd = fusor.critical_domain(
        "preserved",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245", "BRAF", use_minimal_gene_descr=False)
    compare_cd(cd, critical_domain)

    cd = fusor.critical_domain(
        "preserveded",
        "Serine-threonine/tyrosine-protein kinase, catalytic domain",
        "interpro:IPR001245", "BRAF", use_minimal_gene_descr=False)
    assert cd[0] is None
    assert "value is not a valid enumeration member; permitted: " \
           "'lost', 'preserved'" in cd[1]


def test__location_descriptor(fusor, location_descriptor):
    """Test that _location_descriptor method works correctly."""
    ld = fusor._location_descriptor(154170398, 154170399, "NM_152263.3")
    assert ld.dict() == location_descriptor.dict()

    expected = copy.deepcopy(location_descriptor)
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
