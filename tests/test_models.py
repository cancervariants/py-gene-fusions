"""Module for testing the fusion model."""
from pydantic import ValidationError
import pytest
import json
from fusor.models import AbstractFusion, AssayedFusion, CategoricalFusion, \
    TranscriptSegmentComponent, TemplatedSequenceComponent, \
    UnknownGeneComponent, GeneComponent,  AnyGeneComponent, LinkerComponent, \
    FunctionalDomain, Event, RegulatoryElement
import copy
from tests.conftest import EXAMPLES_DIR


@pytest.fixture(scope="module")
def gene_descriptors():
    """Provide possible gene_descriptor input."""
    return [
        {
            "id": "gene:G1",
            "gene": {"gene_id": "hgnc:9339"},
            "label": "G1"
        },
        {
            "id": "gene:ABL",
            "gene": {"gene_id": "hgnc:76"},
            "label": "ABL"
        },
        {
            "id": "gene:BCR1",
            "gene": {"gene_id": "hgnc:1014"},
            "label": "BCR1"
        },
        {
            "id": "gene:NTRK1",
            "gene_id": "hgnc:8031",
            "label": "NTRK1"
        },
        {
            "id": "gene:ALK",
            "gene_id": "hgnc:1837",
            "label": "ALK",
        },
        {
            "id": "gene:YAP1",
            "gene_id": "hgnc:16262",
            "label": "YAP1"
        }
    ]


@pytest.fixture(scope="module")
def location_descriptors():
    """Provide possible templated_sequence input."""
    return [
        {
            "id": "NC_000001.11:15455",
            "type": "LocationDescriptor",
            "location": {
                "sequence_id": "ncbi:NC_000001.11",
                "interval": {
                    "start": {
                        "type": "Number",
                        "value": 15455
                    },
                    "end": {
                        "type": "Number",
                        "value": 15456
                    }
                },
                "type": "SequenceLocation"
            },
            "label": "NC_000001.11:15455",
        },
        {
            "id": "NC_000001.11:15566",
            "type": "LocationDescriptor",
            "location": {
                "sequence_id": "ncbi:NC_000001.11",
                "interval": {
                    "start": {
                        "type": "Number",
                        "value": 15565
                    },
                    "end": {
                        "type": "Number",
                        "value": 15566
                    }
                },
                "type": "SequenceLocation"
            },
            "label": "NC_000001.11:15566",
        },
        {
            "id": "chr12:p12.1",
            "type": "LocationDescriptor",
            "location": {
                "species_id": "taxonomy:9606",
                "chr": "12",
                "interval": {"start": "p12.1", "end": "p12.1"}
            },
            "label": "chr12:p12.1",
        },
        {
            "id": "chr12:p12.2",
            "type": "LocationDescriptor",
            "location": {
                "species_id": "taxonomy:9606",
                "chr": "12",
                "interval": {"start": "p12.2", "end": "p12.2"}
            },
            "label": "chr12:p12.2",
        },
        {
            "id": "NC_000001.11:15455-15566",
            "type": "LocationDescriptor",
            "location": {
                "sequence_id": "ncbi:NC_000001.11",
                "interval": {
                    "start": {
                        "type": "Number",
                        "value": 15455
                    },
                    "end": {
                        "type": "Number",
                        "value": 15566
                    }
                },
                "type": "SequenceLocation"
            },
            "label": "NC_000001.11:15455-15566",
        },
        {
            "id": "chr12:p12.1-p12.2",
            "type": "LocationDescriptor",
            "location": {
                "species_id": "taxonomy:9606",
                "chr": "12",
                "interval": {"start": "p12.1", "end": "p12.2"}
            },
            "label": "chr12:p12.1-p12.2",
        },
        {
            "id": "fusor.location_descriptor:NP_001123617.1",
            "type": "LocationDescriptor",
            "location": {
                "sequence_id": "ga4gh:SQ.sv5egNzqN5koJQH6w0M4tIK9tEDEfJl7",
                "type": "SequenceLocation",
                "interval": {
                    "start": {
                        "type": "Number",
                        "value": 171
                    },
                    "end": {
                        "type": "Number",
                        "value": 204
                    }
                }
            }
        },
        {
            "id": "fusor.location_descriptor:NP_002520.2",
            "type": "LocationDescriptor",
            "location": {
                "sequence_id": "ga4gh:SQ.vJvm06Wl5J7DXHynR9ksW7IK3_3jlFK6",
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
    ]


@pytest.fixture(scope="module")
def functional_domains(gene_descriptors, location_descriptors):
    """Provide possible functional_domains input."""
    return [
        {
            "type": "FunctionalDomain",
            "status": "preserved",
            "name": "WW domain",
            "id": "interpro:IPR001202",
            "gene_descriptor": gene_descriptors[5],
            "location_descriptor": location_descriptors[6]
        },
        {
            "status": "lost",
            "name": "Tyrosine-protein kinase, catalytic domain",
            "id": "interpro:IPR020635",
            "gene_descriptor": gene_descriptors[3],
            "location_descriptor": location_descriptors[7],
        }
    ]


@pytest.fixture(scope="module")
def transcript_segments(location_descriptors, gene_descriptors):
    """Provide possible transcript_segment input."""
    return [
        {
            "transcript": "refseq:NM_152263.3",
            "exon_start": 1,
            "exon_start_offset": -9,
            "exon_end": 8,
            "exon_end_offset": 7,
            "gene_descriptor": gene_descriptors[0],
            "component_genomic_start": location_descriptors[2],
            "component_genomic_end": location_descriptors[3]
        },
        {
            "type": "TranscriptSegmentComponent",
            "transcript": "refseq:NM_034348.3",
            "exon_start": 1,
            "exon_end": 8,
            "gene_descriptor": gene_descriptors[3],
            "component_genomic_start": location_descriptors[0],
            "component_genomic_end": location_descriptors[1]
        },
        {
            "type": "TranscriptSegmentComponent",
            "transcript": "refseq:NM_938439.4",
            "exon_start": 7,
            "exon_end": 14,
            "exon_end_offset": -5,
            "gene_descriptor": gene_descriptors[4],
            "component_genomic_start": location_descriptors[0],
            "component_genomic_end": location_descriptors[1]
        },
        {
            "type": "TranscriptSegmentComponent",
            "transcript": "refseq:NM_938439.4",
            "exon_start": 7,
            "gene_descriptor": gene_descriptors[4],
            "component_genomic_start": location_descriptors[0]
        }
    ]


@pytest.fixture(scope="module")
def gene_components(gene_descriptors):
    """Provide possible gene component input data."""
    return [
        {
            "type": "GeneComponent",
            "gene_descriptor": gene_descriptors[1],
        }
    ]


@pytest.fixture(scope="module")
def templated_sequence_components(location_descriptors):
    """Provide possible templated sequence component input data."""
    return [
        {
            "type": "TemplatedSequenceComponent",
            "strand": "+",
            "region": location_descriptors[5]
        },
        {
            "type": "TemplatedSequenceComponent",
            "strand": "-",
            "region": location_descriptors[4]
        }
    ]


@pytest.fixture(scope="module")
def sequence_descriptors():
    """Provide possible SequenceDescriptor input data"""
    return [
        {
            "id": "sequence:ACGT",
            "type": "SequenceDescriptor",
            "sequence": "ACGT",
            "residue_type": "SO:0000348"
        },
        {
            "id": "sequence:T",
            "type": "SequenceDescriptor",
            "sequence": "T",
            "residue_type": "SO:0000348"
        },
        {
            "id": "sequence:actgu",
            "type": "SequenceDescriptor",
            "sequence": "actgu",
            "residue_type": "SO:0000348"
        }
    ]


@pytest.fixture(scope="module")
def linkers(sequence_descriptors):
    """Provide possible linker component input data."""
    return [
        {
            "type": "LinkerSequenceComponent",
            "linker_sequence": sequence_descriptors[0]
        },
        {
            "type": "LinkerSequenceComponent",
            "linker_sequence": sequence_descriptors[1]
        },
        {
            "type": "LinkerSequenceComponent",
            "linker_sequence": sequence_descriptors[2]
        }
    ]


@pytest.fixture(scope="module")
def regulatory_elements(gene_descriptors):
    """Provide possible regulatory_element input data."""
    return [
        {
            "element_type": "promoter",
            "associated_gene": gene_descriptors[0]
        }
    ]


def check_validation_error(exc_info, expected_msg: str,
                           index: int = 0):
    """Check ValidationError instance for expected message.
    :param ExceptionInfo exc_info: ValidationError instance raised and captured
    by pytest.
    :param str expected_msg: message expected to be provided by error
    :param int index: optional index (if multiple errors are raised)
    :return: None, but may raise AssertionError if incorrect behavior found.
    """
    assert exc_info.value.errors()[index]["msg"] == expected_msg


def test_functional_domain(functional_domains, gene_descriptors):
    """Test FunctionalDomain object initializes correctly"""
    test_domain = FunctionalDomain(**functional_domains[0])
    assert test_domain.type == "FunctionalDomain"
    assert test_domain.status == "preserved"
    assert test_domain.name == "WW domain"
    assert test_domain.id == "interpro:IPR001202"
    assert test_domain.gene_descriptor.id == "gene:YAP1"
    assert test_domain.gene_descriptor.gene_id == "hgnc:16262"
    assert test_domain.gene_descriptor.label == "YAP1"
    test_loc = test_domain.location_descriptor
    assert test_loc.id == "fusor.location_descriptor:NP_001123617.1"
    assert test_loc.type == "LocationDescriptor"
    assert test_loc.location.sequence_id == "ga4gh:SQ.sv5egNzqN5koJQH6w0M4tIK9tEDEfJl7"  # noqa: E501
    assert test_loc.location.interval.type == "SequenceInterval"
    assert test_loc.location.interval.start.value == 171
    assert test_loc.location.interval.end.value == 204

    test_domain = FunctionalDomain(**functional_domains[1])
    assert test_domain.type == "FunctionalDomain"
    assert test_domain.status == "lost"
    assert test_domain.name == "Tyrosine-protein kinase, catalytic domain"
    assert test_domain.id == "interpro:IPR020635"
    assert test_domain.gene_descriptor.id == "gene:NTRK1"
    assert test_domain.gene_descriptor.gene_id == "hgnc:8031"
    assert test_domain.gene_descriptor.label == "NTRK1"
    test_loc = test_domain.location_descriptor
    assert test_loc.id == "fusor.location_descriptor:NP_002520.2"
    assert test_loc.type == "LocationDescriptor"
    assert test_loc.location.sequence_id == "ga4gh:SQ.vJvm06Wl5J7DXHynR9ksW7IK3_3jlFK6"  # noqa: E501
    assert test_loc.location.interval.type == "SequenceInterval"
    assert test_loc.location.interval.start.value == 510
    assert test_loc.location.interval.end.value == 781

    # test status string
    with pytest.raises(ValidationError) as exc_info:
        FunctionalDomain(**{
            "status": "gained",
            "name": "tyrosine kinase catalytic domain",
            "id": "interpro:IPR020635",
            "gene_descriptor": gene_descriptors[0]
        })
    msg = "value is not a valid enumeration member; permitted: 'lost', 'preserved'"  # noqa: E501
    check_validation_error(exc_info, msg)

    # test domain ID CURIE requirement
    with pytest.raises(ValidationError) as exc_info:
        FunctionalDomain(**{
            "status": "lost",
            "name": "tyrosine kinase catalytic domain",
            "id": "interpro_IPR020635",
            "gene_descriptor": gene_descriptors[0]
        })
    msg = "string does not match regex \"^\w[^:]*:.+$\""  # noqa: W605, Q003
    check_validation_error(exc_info, msg)


def test_transcript_segment_component(transcript_segments):
    """Test TranscriptSegmentComponent object initializes correctly"""
    test_component = TranscriptSegmentComponent(**transcript_segments[0])
    assert test_component.transcript == "refseq:NM_152263.3"
    assert test_component.exon_start == 1
    assert test_component.exon_start_offset == -9
    assert test_component.exon_end == 8
    assert test_component.exon_end_offset == 7
    assert test_component.gene_descriptor.id == "gene:G1"
    assert test_component.gene_descriptor.label == "G1"
    assert test_component.gene_descriptor.gene.gene_id == "hgnc:9339"
    test_region_start = test_component.component_genomic_start
    assert test_region_start.location.species_id == "taxonomy:9606"
    assert test_region_start.location.type == "ChromosomeLocation"
    assert test_region_start.location.chr == "12"
    assert test_region_start.location.interval.start == "p12.1"
    assert test_region_start.location.interval.end == "p12.1"
    test_region_end = test_component.component_genomic_end
    assert test_region_end.location.species_id == "taxonomy:9606"
    assert test_region_end.location.type == "ChromosomeLocation"
    assert test_region_end.location.chr == "12"
    assert test_region_end.location.interval.start == "p12.2"
    assert test_region_end.location.interval.end == "p12.2"

    test_component = TranscriptSegmentComponent(**transcript_segments[3])
    assert test_component.transcript == "refseq:NM_938439.4"
    assert test_component.exon_start == 7
    assert test_component.exon_start_offset == 0
    assert test_component.exon_end is None
    assert test_component.exon_end_offset is None

    # check CURIE requirement
    with pytest.raises(ValidationError) as exc_info:
        TranscriptSegmentComponent(**{
            "transcript": "NM_152263.3",
            "exon_start": "1",
            "exon_start_offset": "-9",
            "exon_end": "8",
            "exon_end_offset": "7",
            "gene_descriptor": {
                "id": "test:1",
                "gene": {"id": "hgnc:1"},
                "label": "G1"
            },
            "component_genomic_start": {
                "location": {
                    "species_id": "taxonomy:9606", "chr": "12",
                    "interval": {"start": "p12.1", "end": "p12.1"},
                }
            },
            "component_genomic_end": {
                "location": {
                    "species_id": "taxonomy:9606", "chr": "12",
                    "interval": {"start": "p12.2", "end": "p12.2"},
                }
            }
        })
    msg = "string does not match regex \"^\\w[^:]*:.+$\""  # noqa: W605, Q003
    check_validation_error(exc_info, msg)

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert TranscriptSegmentComponent(**{
            "type": "TemplatedSequenceComponent",
            "transcript": "NM_152263.3",
            "exon_start": "1",
            "exon_start_offset": "-9",
            "exon_end": "8",
            "exon_end_offset": "7",
            "gene_descriptor": {
                "id": "test:1",
                "gene": {"id": "hgnc:1"},
                "label": "G1"
            },
            "component_genomic_start": {
                "location": {
                    "species_id": "taxonomy:9606", "chr": "12",
                    "interval": {"start": "p12.1", "end": "p12.2"},
                }
            },
            "component_genomic_end": {
                "location": {
                    "species_id": "taxonomy:9606", "chr": "12",
                    "interval": {"start": "p12.2", "end": "p12.2"},
                }
            }
        })
    msg = "unexpected value; permitted: <FUSORTypes.TRANSCRIPT_SEGMENT_COMPONENT: 'TranscriptSegmentComponent'>"  # noqa: E501
    check_validation_error(exc_info, msg)

    # test component required
    with pytest.raises(ValidationError) as exc_info:
        assert TranscriptSegmentComponent(**{
            "component_type": "templated_sequence",
            "transcript": "NM_152263.3",
            "exon_start": "1",
            "exon_start_offset": "-9",
            "gene_descriptor": {
                "id": "test:1",
                "gene": {"id": "hgnc:1"},
                "label": "G1"
            }
        })
    msg = "Must give `component_genomic_start` if `exon_start` is given"
    check_validation_error(exc_info, msg)

    # Neither exon_start or exon_end given
    with pytest.raises(ValidationError) as exc_info:
        assert TranscriptSegmentComponent(**{
            "type": "TranscriptSegmentComponent",
            "transcript": "NM_152263.3",
            "exon_start_offset": "-9",
            "exon_end_offset": "7",
            "gene_descriptor": {
                "id": "test:1",
                "gene": {"id": "hgnc:1"},
                "label": "G1"
            },
            "component_genomic_start": {
                "location": {
                    "species_id": "taxonomy:9606", "chr": "12",
                    "interval": {"start": "p12.1", "end": "p12.2"},
                }
            },
            "component_genomic_end": {
                "location": {
                    "species_id": "taxonomy:9606", "chr": "12",
                    "interval": {"start": "p12.2", "end": "p12.2"},
                }
            }
        })
    msg = "Must give values for either `exon_start`, `exon_end`, or both"
    check_validation_error(exc_info, msg)


def test_linker_component(linkers):
    """Test Linker object initializes correctly"""
    def check_linker(actual, expected_id, expected_sequence):
        assert actual.type == "LinkerSequenceComponent"
        assert actual.linker_sequence.id == expected_id
        assert actual.linker_sequence.sequence == expected_sequence
        assert actual.linker_sequence.type == "SequenceDescriptor"
        assert actual.linker_sequence.residue_type == "SO:0000348"

    for args in ((LinkerComponent(**linkers[0]), "sequence:ACGT", "ACGT"),
                 (LinkerComponent(**linkers[1]), "sequence:T", "T"),
                 (LinkerComponent(**linkers[2]), "sequence:actgu", "ACTGU")):
        check_linker(*args)

    # check base validation
    with pytest.raises(ValidationError) as exc_info:
        LinkerComponent(**{
            "linker_sequence":
            {
                "id": "sequence:ACT1",
                "sequence": "ACT1"
            }
        })
    msg = "sequence does not match regex '^[A-Za-z*\\-]*$'"
    check_validation_error(exc_info, msg)

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert LinkerComponent(**{
            "type": "TemplatedSequenceComponent",
            "linker_sequence":
            {
                "id": "sequence:ATG",
                "sequence": "ATG"
            }
        })
    msg = "unexpected value; permitted: <FUSORTypes.LINKER_SEQUENCE_COMPONENT: 'LinkerSequenceComponent'>"  # noqa: E501
    check_validation_error(exc_info, msg)

    # test no extras
    with pytest.raises(ValidationError) as exc_info:
        assert LinkerComponent(**{
            "type": "LinkerSequenceComponent",
            "linker_sequence": {
                "id": "sequence:G",
                "sequence": "G"
            },
            "bonus_value": "bonus",
        })
    msg = "extra fields not permitted"
    check_validation_error(exc_info, msg)


def test_genomic_region_component(templated_sequence_components,
                                  location_descriptors):
    """Test that TemplatedSequenceComponent initializes correctly."""
    def assert_genomic_region_test_component(test):
        """Assert that test templated_sequence_components[0] data matches
        expected values.
        """
        assert test.type == "TemplatedSequenceComponent"
        assert test.strand.value == "+"
        assert test.region.id == "chr12:p12.1-p12.2"
        assert test.region.type == "LocationDescriptor"
        assert test.region.location.species_id == "taxonomy:9606"
        assert test.region.location.chr == "12"
        assert test.region.location.interval.start == "p12.1"
        assert test.region.location.interval.end == "p12.2"
        assert test.region.label == "chr12:p12.1-p12.2"

    test_component = \
        TemplatedSequenceComponent(**templated_sequence_components[0])
    assert_genomic_region_test_component(test_component)

    genomic_region_components_cpy = \
        copy.deepcopy(templated_sequence_components[0])
    genomic_region_components_cpy["region"]["location"]["_id"] = "location:1"
    test_component = \
        TemplatedSequenceComponent(**genomic_region_components_cpy)
    assert_genomic_region_test_component(test_component)

    genomic_region_components_cpy = \
        copy.deepcopy(templated_sequence_components[0])
    genomic_region_components_cpy["region"]["location_id"] = "location:1"
    test_component = \
        TemplatedSequenceComponent(**genomic_region_components_cpy)
    assert_genomic_region_test_component(test_component)

    with pytest.raises(ValidationError) as exc_info:
        TemplatedSequenceComponent(**{
            "region": {
                "interval": {
                    "start": 39408,
                    "stop": 39414
                }
            },
            "sequence_id": "ga4gh:SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl"
        })
    msg = "Must give values for either `location`, `location_id`, or both"
    check_validation_error(exc_info, msg)

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert TemplatedSequenceComponent(**{
            "type": "GeneComponent",
            "region": location_descriptors[0],
            "strand": "+"
        })
    msg = "unexpected value; permitted: <FUSORTypes.TEMPLATED_SEQUENCE_COMPONENT: 'TemplatedSequenceComponent'>"  # noqa: E501
    check_validation_error(exc_info, msg)


def test_gene_component(gene_descriptors):
    """Test that Gene component initializes correctly."""
    test_component = GeneComponent(**{"gene_descriptor": gene_descriptors[0]})
    assert test_component.type == "GeneComponent"
    assert test_component.gene_descriptor.id == "gene:G1"
    assert test_component.gene_descriptor.label == "G1"
    assert test_component.gene_descriptor.gene.gene_id == "hgnc:9339"

    # test CURIE requirement
    with pytest.raises(ValidationError) as exc_info:
        GeneComponent(**{
            "gene_descriptor": {
                "id": "G1",
                "gene": {"gene_id": "hgnc:9339"},
                "label": "G1"
            }
        })
    msg = "string does not match regex \"^\w[^:]*:.+$\""  # noqa: W605, Q003
    check_validation_error(exc_info, msg)

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert GeneComponent(**{
            "type": "UnknownGeneComponent",
            "gene_descriptor": gene_descriptors[0]
        })
    msg = "unexpected value; permitted: <FUSORTypes.GENE_COMPONENT: 'GeneComponent'>"  # noqa: E501
    check_validation_error(exc_info, msg)


def test_unknown_gene_component():
    """Test that unknown_gene component initializes correctly."""
    test_component = UnknownGeneComponent()
    assert test_component.type == "UnknownGeneComponent"

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert UnknownGeneComponent(type="gene")
    msg = "unexpected value; permitted: <FUSORTypes.UNKNOWN_GENE_COMPONENT: 'UnknownGeneComponent'>"  # noqa: E501
    check_validation_error(exc_info, msg)


def test_any_gene_component():
    """Test that any_gene component initializes correctly."""
    test_component = AnyGeneComponent()
    assert test_component.type == "AnyGeneComponent"

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert AnyGeneComponent(type="unknown_gene")
    msg = "unexpected value; permitted: <FUSORTypes.ANY_GENE_COMPONENT: 'AnyGeneComponent'>"  # noqa: E501
    check_validation_error(exc_info, msg)


def test_event():
    """Test Event object initializes correctly"""
    rearrangement = "rearrangement"
    test_event = Event(rearrangement)
    assert test_event.value == rearrangement

    with pytest.raises(ValueError):
        Event("combination")


def test_regulatory_element(regulatory_elements, gene_descriptors):
    """Test RegulatoryElement object initializes correctly"""
    test_reg_elmt = RegulatoryElement(**regulatory_elements[0])
    assert test_reg_elmt.element_type.value == "promoter"
    assert test_reg_elmt.associated_gene.id == "gene:G1"
    assert test_reg_elmt.associated_gene.gene.gene_id == "hgnc:9339"
    assert test_reg_elmt.associated_gene.label == "G1"

    # check type constraint
    with pytest.raises(ValidationError) as exc_info:
        RegulatoryElement(**{
            "element_type": "notpromoter",
            "associated_gene": gene_descriptors[0]
        })
    assert exc_info.value.errors()[0]["msg"].startswith(
        "value is not a valid enumeration member; permitted: "
    )

    # require minimum input
    with pytest.raises(ValidationError) as exc_info:
        RegulatoryElement(**{
            "element_type": "enhancer",
        })
    assert exc_info.value.errors()[0]["msg"] == "Must set >=1 of {`element_reference`, `associated_gene`, `genomic_location`}"  # noqa: E501


def test_fusion(functional_domains, transcript_segments,
                templated_sequence_components, linkers, gene_components,
                regulatory_elements):
    """Test that Fusion object initializes correctly"""
    unknown_component = {
        "type": "UnknownGeneComponent",
    }

    # test valid object
    fusion = CategoricalFusion(**{
        "r_frame_preserved": True,
        "functional_domains": [functional_domains[0]],
        "structural_components": [
            transcript_segments[1], transcript_segments[2]
        ],
        "regulatory_elements": [regulatory_elements[0]]
    })

    assert fusion.structural_components[0].transcript == "refseq:NM_034348.3"

    # check correct parsing of nested items
    fusion = CategoricalFusion(**{
        "structural_components": [
            {
                "type": "GeneComponent",
                "gene_descriptor": {
                    "type": "GeneDescriptor",
                    "id": "gene:NTRK1",
                    "label": "NTRK1",
                    "gene_id": "hgnc:8031"
                }
            },
            {
                "type": "GeneComponent",
                "gene_descriptor": {
                    "type": "GeneDescriptor",
                    "id": "gene:ABL1",
                    "label": "ABL1",
                    "gene_id": "hgnc:76"
                }
            }
        ],
        "regulatory_elements": []
    })
    assert fusion.structural_components[0].type == "GeneComponent"
    assert fusion.structural_components[0].gene_descriptor.id == "gene:NTRK1"
    assert fusion.structural_components[1].type == "GeneComponent"
    assert fusion.structural_components[1].gene_descriptor.type == "GeneDescriptor"  # noqa: E501

    # test that non-component properties are optional
    assert CategoricalFusion(**{
        "structural_components": [
            transcript_segments[1], transcript_segments[2]
        ]
    })

    # test variety of component types
    assert AssayedFusion(**{
        "type": "AssayedFusion",
        "structural_components": [
            unknown_component,
            gene_components[0],
            transcript_segments[2],
            templated_sequence_components[1],
            linkers[0],
        ]
    })
    with pytest.raises(ValidationError) as exc_info:
        assert CategoricalFusion(**{
            "type": "CategoricalFusion",
            "structural_components": [
                {
                    "type": "LinkerSequenceComponent",
                    "linker_sequence": {
                        "id": "a:b",
                        "type": "SequenceDescriptor",
                        "sequence": "AC",
                        "residue_type": "SO:0000348"
                    }
                },
                {
                    "type": "LinkerSequenceComponent",
                    "linker_sequence": {
                        "id": "a:b",
                        "type": "SequenceDescriptor",
                        "sequence": "AC",
                        "residue_type": "SO:0000348"
                    }
                }
            ]
        })
    msg = "First structural component cannot be LinkerSequence"
    check_validation_error(exc_info, msg)

    # components are mandatory
    with pytest.raises(ValidationError) as exc_info:
        assert AssayedFusion(**{
            "functional_domains": [functional_domains[1]],
            "causative_event": "rearrangement",
            "regulatory_elements": [regulatory_elements[0]]
        })
    check_validation_error(exc_info, "field required")

    # must have >= 2 components + regulatory elements
    with pytest.raises(ValidationError) as exc_info:
        assert AssayedFusion(**{
            "structural_components": [unknown_component]
        })
    msg = "Provided fusion contains an insufficient number of structural components and regulatory elements."  # noqa: E501
    check_validation_error(exc_info, msg)
    assert AssayedFusion(**{
        "regulatory_elements": [regulatory_elements[0]],
        "structural_components": [transcript_segments[0]]
    })

    # can't create base fusion
    with pytest.raises(ValidationError) as exc_info:
        assert AbstractFusion(**{
            "structural_components": [transcript_segments[2],
                                      linkers[0]]
        })
    check_validation_error(exc_info,
                           "Cannot instantiate Fusion abstract class")


def test_examples(exhaustive_example, fusion_example):
    """Test example JSON files."""
    assert CategoricalFusion(**exhaustive_example)

    assert CategoricalFusion(**fusion_example)

    with open(EXAMPLES_DIR / "alk.json", "r") as alk_file:
        assert CategoricalFusion(**json.load(alk_file))

    with open(EXAMPLES_DIR / "epcam_msh2.json", "r") as epcam_file:
        assert CategoricalFusion(**json.load(epcam_file))

    with open(EXAMPLES_DIR / "tpm3_ntrk1.json", "r") as ntrk_file:
        assert AssayedFusion(**json.load(ntrk_file))

    with open(EXAMPLES_DIR / "tpm3_pdgfrb.json", "r") as pdgfrb_file:
        assert CategoricalFusion(**json.load(pdgfrb_file))

    with open(EXAMPLES_DIR / "ewsr1.json", "r") as ewsr1_file:
        assert AssayedFusion(**json.load(ewsr1_file))
