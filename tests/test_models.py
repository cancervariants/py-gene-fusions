"""Module for testing the fusion model."""
import copy

import pytest
from pydantic import ValidationError

from fusor.models import (
    AbstractFusion,
    Assay,
    AssayedFusion,
    CategoricalFusion,
    CausativeEvent,
    EventType,
    FunctionalDomain,
    GeneElement,
    LinkerElement,
    MultiplePossibleGenesElement,
    RegulatoryElement,
    TemplatedSequenceElement,
    TranscriptSegmentElement,
    UnknownGeneElement,
)


@pytest.fixture(scope="module")
def gene_descriptors():
    """Provide possible gene_descriptor input."""
    return [
        {"id": "gene:G1", "gene": {"gene_id": "hgnc:9339"}, "label": "G1"},
        {"id": "gene:ABL", "gene": {"gene_id": "hgnc:76"}, "label": "ABL"},
        {"id": "gene:BCR1", "gene": {"gene_id": "hgnc:1014"}, "label": "BCR1"},
        {"id": "gene:NTRK1", "gene_id": "hgnc:8031", "label": "NTRK1"},
        {
            "id": "gene:ALK",
            "gene_id": "hgnc:1837",
            "label": "ALK",
        },
        {"id": "gene:YAP1", "gene_id": "hgnc:16262", "label": "YAP1"},
        # alternate structure
        {
            "id": "normalize.gene:BRAF",
            "type": "GeneDescriptor",
            "label": "BRAF",
            "gene_id": "hgnc:1097",
        },
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
                    "start": {"type": "Number", "value": 15455},
                    "end": {"type": "Number", "value": 15456},
                },
                "type": "SequenceLocation",
            },
            "label": "NC_000001.11:15455",
        },
        {
            "id": "NC_000001.11:15566",
            "type": "LocationDescriptor",
            "location": {
                "sequence_id": "ncbi:NC_000001.11",
                "interval": {
                    "start": {"type": "Number", "value": 15565},
                    "end": {"type": "Number", "value": 15566},
                },
                "type": "SequenceLocation",
            },
            "label": "NC_000001.11:15566",
        },
        {
            "id": "chr12:p12.1",
            "type": "LocationDescriptor",
            "location": {
                "species_id": "taxonomy:9606",
                "chr": "12",
                "interval": {"start": "p12.1", "end": "p12.1"},
            },
            "label": "chr12:p12.1",
        },
        {
            "id": "chr12:p12.2",
            "type": "LocationDescriptor",
            "location": {
                "species_id": "taxonomy:9606",
                "chr": "12",
                "interval": {"start": "p12.2", "end": "p12.2"},
            },
            "label": "chr12:p12.2",
        },
        {
            "id": "NC_000001.11:15455-15566",
            "type": "LocationDescriptor",
            "location": {
                "sequence_id": "ncbi:NC_000001.11",
                "interval": {
                    "start": {"type": "Number", "value": 15455},
                    "end": {"type": "Number", "value": 15566},
                },
                "type": "SequenceLocation",
            },
            "label": "NC_000001.11:15455-15566",
        },
        {
            "id": "chr12:p12.1-p12.2",
            "type": "LocationDescriptor",
            "location": {
                "species_id": "taxonomy:9606",
                "chr": "12",
                "interval": {"start": "p12.1", "end": "p12.2"},
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
                    "start": {"type": "Number", "value": 171},
                    "end": {"type": "Number", "value": 204},
                },
            },
        },
        {
            "id": "fusor.location_descriptor:NP_002520.2",
            "type": "LocationDescriptor",
            "location": {
                "sequence_id": "ga4gh:SQ.vJvm06Wl5J7DXHynR9ksW7IK3_3jlFK6",
                "type": "SequenceLocation",
                "interval": {
                    "start": {"type": "Number", "value": 510},
                    "end": {"type": "Number", "value": 781},
                },
            },
        },
    ]


@pytest.fixture(scope="module")
def functional_domains(gene_descriptors, location_descriptors):
    """Provide possible functional_domains input."""
    return [
        {
            "type": "FunctionalDomain",
            "status": "preserved",
            "label": "WW domain",
            "_id": "interpro:IPR001202",
            "associated_gene": gene_descriptors[5],
            "sequence_location": location_descriptors[6],
        },
        {
            "status": "lost",
            "label": "Tyrosine-protein kinase, catalytic domain",
            "_id": "interpro:IPR020635",
            "associated_gene": gene_descriptors[3],
            "sequence_location": location_descriptors[7],
        },
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
            "element_genomic_start": location_descriptors[2],
            "element_genomic_end": location_descriptors[3],
        },
        {
            "type": "TranscriptSegmentElement",
            "transcript": "refseq:NM_034348.3",
            "exon_start": 1,
            "exon_end": 8,
            "gene_descriptor": gene_descriptors[3],
            "element_genomic_start": location_descriptors[0],
            "element_genomic_end": location_descriptors[1],
        },
        {
            "type": "TranscriptSegmentElement",
            "transcript": "refseq:NM_938439.4",
            "exon_start": 7,
            "exon_end": 14,
            "exon_end_offset": -5,
            "gene_descriptor": gene_descriptors[4],
            "element_genomic_start": location_descriptors[0],
            "element_genomic_end": location_descriptors[1],
        },
        {
            "type": "TranscriptSegmentElement",
            "transcript": "refseq:NM_938439.4",
            "exon_start": 7,
            "gene_descriptor": gene_descriptors[4],
            "element_genomic_start": location_descriptors[0],
        },
    ]


@pytest.fixture(scope="module")
def gene_elements(gene_descriptors):
    """Provide possible gene element input data."""
    return [
        {
            "type": "GeneElement",
            "gene_descriptor": gene_descriptors[1],
        },
        {"type": "GeneElement", "gene_descriptor": gene_descriptors[0]},
        {"type"},
    ]


@pytest.fixture(scope="module")
def templated_sequence_elements(location_descriptors):
    """Provide possible templated sequence element input data."""
    return [
        {
            "type": "TemplatedSequenceElement",
            "strand": "+",
            "region": location_descriptors[5],
        },
        {
            "type": "TemplatedSequenceElement",
            "strand": "-",
            "region": location_descriptors[4],
        },
    ]


@pytest.fixture(scope="module")
def sequence_descriptors():
    """Provide possible SequenceDescriptor input data"""
    return [
        {
            "id": "sequence:ACGT",
            "type": "SequenceDescriptor",
            "sequence": "ACGT",
            "residue_type": "SO:0000348",
        },
        {
            "id": "sequence:T",
            "type": "SequenceDescriptor",
            "sequence": "T",
            "residue_type": "SO:0000348",
        },
        {
            "id": "sequence:actgu",
            "type": "SequenceDescriptor",
            "sequence": "actgu",
            "residue_type": "SO:0000348",
        },
    ]


@pytest.fixture(scope="module")
def linkers(sequence_descriptors):
    """Provide possible linker element input data."""
    return [
        {"type": "LinkerSequenceElement", "linker_sequence": sequence_descriptors[0]},
        {"type": "LinkerSequenceElement", "linker_sequence": sequence_descriptors[1]},
        {"type": "LinkerSequenceElement", "linker_sequence": sequence_descriptors[2]},
    ]


@pytest.fixture(scope="module")
def unknown_element():
    """Provide UnknownGene element."""
    return {"type": "UnknownGeneElement"}


@pytest.fixture(scope="module")
def regulatory_elements(gene_descriptors):
    """Provide possible regulatory_element input data."""
    return [{"regulatory_class": "promoter", "associated_gene": gene_descriptors[0]}]


def check_validation_error(exc_info, expected_msg: str, index: int = 0):
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
    assert test_domain.label == "WW domain"
    assert test_domain.id == "interpro:IPR001202"
    assert test_domain.associated_gene.id == "gene:YAP1"
    assert test_domain.associated_gene.gene_id == "hgnc:16262"
    assert test_domain.associated_gene.label == "YAP1"
    test_loc = test_domain.sequence_location
    assert test_loc.id == "fusor.location_descriptor:NP_001123617.1"
    assert test_loc.type == "LocationDescriptor"
    assert test_loc.location.sequence_id == "ga4gh:SQ.sv5egNzqN5koJQH6w0M4tIK9tEDEfJl7"
    assert test_loc.location.interval.type == "SequenceInterval"
    assert test_loc.location.interval.start.value == 171
    assert test_loc.location.interval.end.value == 204

    test_domain = FunctionalDomain(**functional_domains[1])
    assert test_domain.type == "FunctionalDomain"
    assert test_domain.status == "lost"
    assert test_domain.label == "Tyrosine-protein kinase, catalytic domain"
    assert test_domain.id == "interpro:IPR020635"
    assert test_domain.associated_gene.id == "gene:NTRK1"
    assert test_domain.associated_gene.gene_id == "hgnc:8031"
    assert test_domain.associated_gene.label == "NTRK1"
    test_loc = test_domain.sequence_location
    assert test_loc.id == "fusor.location_descriptor:NP_002520.2"
    assert test_loc.type == "LocationDescriptor"
    assert test_loc.location.sequence_id == "ga4gh:SQ.vJvm06Wl5J7DXHynR9ksW7IK3_3jlFK6"
    assert test_loc.location.interval.type == "SequenceInterval"
    assert test_loc.location.interval.start.value == 510
    assert test_loc.location.interval.end.value == 781

    # test status string
    with pytest.raises(ValidationError) as exc_info:
        FunctionalDomain(
            status="gained",
            name="tyrosine kinase catalytic domain",
            id="interpro:IPR020635",
            associated_gene=gene_descriptors[0],
        )
    msg = "Input should be 'lost' or 'preserved'"
    check_validation_error(exc_info, msg)

    # test domain ID CURIE requirement
    with pytest.raises(ValidationError) as exc_info:
        FunctionalDomain(
            status="lost",
            label="tyrosine kinase catalytic domain",
            id="interpro_IPR020635",
            associated_gene=gene_descriptors[0],
        )
    msg = "String should match pattern '^\\w[^:]*:.+$'"
    check_validation_error(exc_info, msg)


def test_transcript_segment_element(transcript_segments):
    """Test TranscriptSegmentElement object initializes correctly"""
    test_element = TranscriptSegmentElement(**transcript_segments[0])
    assert test_element.transcript == "refseq:NM_152263.3"
    assert test_element.exon_start == 1
    assert test_element.exon_start_offset == -9
    assert test_element.exon_end == 8
    assert test_element.exon_end_offset == 7
    assert test_element.gene_descriptor.id == "gene:G1"
    assert test_element.gene_descriptor.label == "G1"
    assert test_element.gene_descriptor.gene.gene_id == "hgnc:9339"
    test_region_start = test_element.element_genomic_start
    assert test_region_start.location.species_id == "taxonomy:9606"
    assert test_region_start.location.type == "ChromosomeLocation"
    assert test_region_start.location.chr == "12"
    assert test_region_start.location.interval.start == "p12.1"
    assert test_region_start.location.interval.end == "p12.1"
    test_region_end = test_element.element_genomic_end
    assert test_region_end.location.species_id == "taxonomy:9606"
    assert test_region_end.location.type == "ChromosomeLocation"
    assert test_region_end.location.chr == "12"
    assert test_region_end.location.interval.start == "p12.2"
    assert test_region_end.location.interval.end == "p12.2"

    test_element = TranscriptSegmentElement(**transcript_segments[3])
    assert test_element.transcript == "refseq:NM_938439.4"
    assert test_element.exon_start == 7
    assert test_element.exon_start_offset == 0
    assert test_element.exon_end is None
    assert test_element.exon_end_offset is None

    # check CURIE requirement
    with pytest.raises(ValidationError) as exc_info:
        TranscriptSegmentElement(
            transcript="NM_152263.3",
            exon_start="1",
            exon_start_offset="-9",
            exon_end="8",
            exon_end_offset="7",
            gene_descriptor={
                "id": "test:1",
                "gene": {"id": "hgnc:1"},
                "label": "G1",
            },
            element_genomic_start={
                "location": {
                    "species_id": "taxonomy:9606",
                    "chr": "12",
                    "interval": {"start": "p12.1", "end": "p12.1"},
                }
            },
            element_genomic_end={
                "location": {
                    "species_id": "taxonomy:9606",
                    "chr": "12",
                    "interval": {"start": "p12.2", "end": "p12.2"},
                }
            },
        )
    msg = "String should match pattern '^\\w[^:]*:.+$'"
    check_validation_error(exc_info, msg)

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert TranscriptSegmentElement(
            type="TemplatedSequenceElement",
            transcript="NM_152263.3",
            exon_start="1",
            exon_start_offset="-9",
            exon_end="8",
            exon_end_offset="7",
            gene_descriptor={
                "id": "test:1",
                "gene": {"id": "hgnc:1"},
                "label": "G1",
            },
            element_genomic_start={
                "location": {
                    "species_id": "taxonomy:9606",
                    "chr": "12",
                    "interval": {"start": "p12.1", "end": "p12.2"},
                }
            },
            element_genomic_end={
                "location": {
                    "species_id": "taxonomy:9606",
                    "chr": "12",
                    "interval": {"start": "p12.2", "end": "p12.2"},
                }
            },
        )
    msg = "Input should be <FUSORTypes.TRANSCRIPT_SEGMENT_ELEMENT: 'TranscriptSegmentElement'>"
    check_validation_error(exc_info, msg)

    # test element required
    with pytest.raises(ValidationError) as exc_info:
        assert TranscriptSegmentElement(
            element_type="templated_sequence",
            transcript="NM_152263.3",
            exon_start="1",
            exon_start_offset="-9",
            gene_descriptor={
                "id": "test:1",
                "gene": {"id": "hgnc:1"},
                "label": "G1",
            },
        )
    msg = "Value error, Must give `element_genomic_start` if `exon_start` is given"
    check_validation_error(exc_info, msg)

    # Neither exon_start or exon_end given
    with pytest.raises(ValidationError) as exc_info:
        assert TranscriptSegmentElement(
            type="TranscriptSegmentElement",
            transcript="NM_152263.3",
            exon_start_offset="-9",
            exon_end_offset="7",
            gene_descriptor={
                "id": "test:1",
                "gene": {"id": "hgnc:1"},
                "label": "G1",
            },
            element_genomic_start={
                "location": {
                    "species_id": "taxonomy:9606",
                    "chr": "12",
                    "interval": {"start": "p12.1", "end": "p12.2"},
                }
            },
            element_genomic_end={
                "location": {
                    "species_id": "taxonomy:9606",
                    "chr": "12",
                    "interval": {"start": "p12.2", "end": "p12.2"},
                }
            },
        )
    msg = "Value error, Must give values for either `exon_start`, `exon_end`, or both"
    check_validation_error(exc_info, msg)


def test_linker_element(linkers):
    """Test Linker object initializes correctly"""

    def check_linker(actual, expected_id, expected_sequence):
        assert actual.type == "LinkerSequenceElement"
        assert actual.linker_sequence.id == expected_id
        assert actual.linker_sequence.sequence == expected_sequence
        assert actual.linker_sequence.type == "SequenceDescriptor"
        assert actual.linker_sequence.residue_type == "SO:0000348"

    for args in (
        (LinkerElement(**linkers[0]), "sequence:ACGT", "ACGT"),
        (LinkerElement(**linkers[1]), "sequence:T", "T"),
        (LinkerElement(**linkers[2]), "sequence:actgu", "ACTGU"),
    ):
        check_linker(*args)

    # check base validation
    with pytest.raises(ValidationError) as exc_info:
        LinkerElement(linker_sequence={"id": "sequence:ACT1", "sequence": "ACT1"})
    msg = "String should match pattern '^[A-Z*\\-]*$'"
    check_validation_error(exc_info, msg)

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert LinkerElement(
            type="TemplatedSequenceElement",
            linker_sequence={"id": "sequence:ATG", "sequence": "ATG"},
        )
    msg = (
        "Input should be <FUSORTypes.LINKER_SEQUENCE_ELEMENT: 'LinkerSequenceElement'>"
    )
    check_validation_error(exc_info, msg)

    # test no extras
    with pytest.raises(ValidationError) as exc_info:
        assert LinkerElement(
            type="LinkerSequenceElement",
            linker_sequence={"id": "sequence:G", "sequence": "G"},
            bonus_value="bonus",
        )
    msg = "Extra inputs are not permitted"
    check_validation_error(exc_info, msg)


def test_genomic_region_element(templated_sequence_elements, location_descriptors):
    """Test that TemplatedSequenceElement initializes correctly."""

    def assert_genomic_region_test_element(test):
        """Assert that test templated_sequence_elements[0] data matches
        expected values.
        """
        assert test.type == "TemplatedSequenceElement"
        assert test.strand.value == "+"
        assert test.region.id == "chr12:p12.1-p12.2"
        assert test.region.type == "LocationDescriptor"
        assert test.region.location.species_id == "taxonomy:9606"
        assert test.region.location.chr == "12"
        assert test.region.location.interval.start == "p12.1"
        assert test.region.location.interval.end == "p12.2"
        assert test.region.label == "chr12:p12.1-p12.2"

    test_element = TemplatedSequenceElement(**templated_sequence_elements[0])
    assert_genomic_region_test_element(test_element)

    genomic_region_elements_cpy = copy.deepcopy(templated_sequence_elements[0])
    genomic_region_elements_cpy["region"]["location"]["_id"] = "location:1"
    test_element = TemplatedSequenceElement(**genomic_region_elements_cpy)
    assert_genomic_region_test_element(test_element)

    genomic_region_elements_cpy = copy.deepcopy(templated_sequence_elements[0])
    genomic_region_elements_cpy["region"]["location_id"] = "location:1"
    test_element = TemplatedSequenceElement(**genomic_region_elements_cpy)
    assert_genomic_region_test_element(test_element)

    with pytest.raises(ValidationError) as exc_info:
        TemplatedSequenceElement(
            region={"interval": {"start": 39408, "stop": 39414}},
            sequence_id="ga4gh:SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl",
        )
    msg = "Field required"
    check_validation_error(exc_info, msg)

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert TemplatedSequenceElement(
            type="GeneElement", region=location_descriptors[0], strand="+"
        )
    msg = "Input should be <FUSORTypes.TEMPLATED_SEQUENCE_ELEMENT: 'TemplatedSequenceElement'>"
    check_validation_error(exc_info, msg)


def test_gene_element(gene_descriptors):
    """Test that Gene Element initializes correctly."""
    test_element = GeneElement(gene_descriptor=gene_descriptors[0])
    assert test_element.type == "GeneElement"
    assert test_element.gene_descriptor.id == "gene:G1"
    assert test_element.gene_descriptor.label == "G1"
    assert test_element.gene_descriptor.gene.gene_id == "hgnc:9339"

    # test CURIE requirement
    with pytest.raises(ValidationError) as exc_info:
        GeneElement(
            gene_descriptor={
                "id": "G1",
                "gene": {"gene_id": "hgnc:9339"},
                "label": "G1",
            }
        )
    msg = "String should match pattern '^\\w[^:]*:.+$'"
    check_validation_error(exc_info, msg)

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert GeneElement(
            type="UnknownGeneElement", gene_descriptor=gene_descriptors[0]
        )
    msg = "Input should be <FUSORTypes.GENE_ELEMENT: 'GeneElement'>"
    check_validation_error(exc_info, msg)


def test_unknown_gene_element():
    """Test that unknown_gene element initializes correctly."""
    test_element = UnknownGeneElement()
    assert test_element.type == "UnknownGeneElement"

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert UnknownGeneElement(type="gene")
    msg = "Input should be <FUSORTypes.UNKNOWN_GENE_ELEMENT: 'UnknownGeneElement'>"
    check_validation_error(exc_info, msg)


def test_mult_gene_element():
    """Test that mult_gene_element initializes correctly."""
    test_element = MultiplePossibleGenesElement()
    assert test_element.type == "MultiplePossibleGenesElement"

    # test enum validation
    with pytest.raises(ValidationError) as exc_info:
        assert MultiplePossibleGenesElement(type="unknown_gene")
    msg = "Input should be <FUSORTypes.MULTIPLE_POSSIBLE_GENES_ELEMENT: 'MultiplePossibleGenesElement'>"
    check_validation_error(exc_info, msg)


def test_event():
    """Test Event object initializes correctly"""
    rearrangement = EventType.REARRANGEMENT
    test_event = CausativeEvent(event_type=rearrangement, event_description=None)
    assert test_event.event_type == rearrangement

    with pytest.raises(ValueError):  # noqa: PT011
        CausativeEvent(event_type="combination")


def test_regulatory_element(regulatory_elements, gene_descriptors):
    """Test RegulatoryElement object initializes correctly"""
    test_reg_elmt = RegulatoryElement(**regulatory_elements[0])
    assert test_reg_elmt.regulatory_class.value == "promoter"
    assert test_reg_elmt.associated_gene.id == "gene:G1"
    assert test_reg_elmt.associated_gene.gene.gene_id == "hgnc:9339"
    assert test_reg_elmt.associated_gene.label == "G1"

    # check type constraint
    with pytest.raises(ValidationError) as exc_info:
        RegulatoryElement(
            regulatory_class="notpromoter", associated_gene=gene_descriptors[0]
        )
    assert exc_info.value.errors()[0]["msg"].startswith("Input should be")

    # require minimum input
    with pytest.raises(ValidationError) as exc_info:
        RegulatoryElement(regulatory_class="enhancer")
    assert (
        exc_info.value.errors()[0]["msg"]
        == "Value error, Must set 1 of {`feature_id`, `associated_gene`} and/or `feature_location`"
    )


def test_fusion(
    functional_domains,
    transcript_segments,
    templated_sequence_elements,
    linkers,
    gene_elements,
    regulatory_elements,
    unknown_element,
):
    """Test that Fusion object initializes correctly"""
    # test valid object
    fusion = CategoricalFusion(
        r_frame_preserved=True,
        critical_functional_domains=[functional_domains[0]],
        structural_elements=[transcript_segments[1], transcript_segments[2]],
        regulatory_element=regulatory_elements[0],
    )

    assert fusion.structural_elements[0].transcript == "refseq:NM_034348.3"

    # check correct parsing of nested items
    fusion = CategoricalFusion(
        structural_elements=[
            {
                "type": "GeneElement",
                "gene_descriptor": {
                    "type": "GeneDescriptor",
                    "id": "gene:NTRK1",
                    "label": "NTRK1",
                    "gene_id": "hgnc:8031",
                },
            },
            {
                "type": "GeneElement",
                "gene_descriptor": {
                    "type": "GeneDescriptor",
                    "id": "gene:ABL1",
                    "label": "ABL1",
                    "gene_id": "hgnc:76",
                },
            },
        ],
        regulatory_element=None,
    )
    assert fusion.structural_elements[0].type == "GeneElement"
    assert fusion.structural_elements[0].gene_descriptor.id == "gene:NTRK1"
    assert fusion.structural_elements[1].type == "GeneElement"
    assert fusion.structural_elements[1].gene_descriptor.type == "GeneDescriptor"

    # test that non-element properties are optional
    assert CategoricalFusion(
        structural_elements=[transcript_segments[1], transcript_segments[2]]
    )

    # test variety of element types
    assert AssayedFusion(
        type="AssayedFusion",
        structural_elements=[
            unknown_element,
            gene_elements[0],
            transcript_segments[2],
            templated_sequence_elements[1],
            linkers[0],
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
    with pytest.raises(ValidationError) as exc_info:
        assert CategoricalFusion(
            type="CategoricalFusion",
            structural_elements=[
                {
                    "type": "LinkerSequenceElement",
                    "linker_sequence": {
                        "id": "a:b",
                        "type": "SequenceDescriptor",
                        "sequence": "AC",
                        "residue_type": "SO:0000348",
                    },
                },
                {
                    "type": "LinkerSequenceElement",
                    "linker_sequence": {
                        "id": "a:b",
                        "type": "SequenceDescriptor",
                        "sequence": "AC",
                        "residue_type": "SO:0000348",
                    },
                },
            ],
        )
    msg = "Value error, First structural element cannot be LinkerSequence"
    check_validation_error(exc_info, msg)


def test_fusion_element_count(
    functional_domains,
    regulatory_elements,
    unknown_element,
    gene_elements,
    transcript_segments,
    gene_descriptors,
):
    """Test fusion element count requirements."""
    # elements are mandatory
    with pytest.raises(ValidationError) as exc_info:
        assert AssayedFusion(
            functional_domains=[functional_domains[1]],
            causative_event="rearrangement",
            regulatory_elements=[regulatory_elements[0]],
        )
    element_ct_msg = (
        "Value error, Fusions must contain >= 2 structural elements, or >=1 structural element "
        "and a regulatory element"
    )
    check_validation_error(exc_info, element_ct_msg)

    # must have >= 2 elements + regulatory elements
    with pytest.raises(ValidationError) as exc_info:
        assert AssayedFusion(
            structural_elements=[unknown_element],
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
    check_validation_error(exc_info, element_ct_msg)

    # unique gene requirements
    uq_gene_error_msg = "Value error, Fusions must form a chimeric transcript from two or more genes, or a novel interaction between a rearranged regulatory element with the expressed product of a partner gene."
    with pytest.raises(ValidationError) as exc_info:
        assert CategoricalFusion(
            structural_elements=[gene_elements[0], gene_elements[0]]
        )
    check_validation_error(exc_info, uq_gene_error_msg)

    with pytest.raises(ValidationError) as exc_info:
        assert CategoricalFusion(
            structural_elements=[gene_elements[1], transcript_segments[0]]
        )
    check_validation_error(exc_info, uq_gene_error_msg)

    with pytest.raises(ValidationError) as exc_info:
        assert CategoricalFusion(
            regulatory_element=regulatory_elements[0],
            structural_elements=[transcript_segments[0]],
        )
    check_validation_error(exc_info, uq_gene_error_msg)

    # use alternate gene descriptor structure
    with pytest.raises(ValidationError) as exc_info:
        assert AssayedFusion(
            type="AssayedFusion",
            structural_elements=[
                {"type": "GeneElement", "gene_descriptor": gene_descriptors[6]},
                {"type": "GeneElement", "gene_descriptor": gene_descriptors[6]},
            ],
            causative_event={
                "type": "CausativeEvent",
                "event_type": "read-through",
            },
            assay={
                "type": "Assay",
                "method_uri": "pmid:33576979",
                "assay_id": "obi:OBI_0003094",
                "assay_name": "fluorescence in-situ hybridization assay",
                "fusion_detection": "inferred",
            },
        )
    with pytest.raises(ValidationError) as exc_info:
        assert AssayedFusion(
            type="AssayedFusion",
            structural_elements=[
                {"type": "GeneElement", "gene_descriptor": gene_descriptors[6]},
            ],
            regulatory_element={
                "type": "RegulatoryElement",
                "regulatory_class": "enhancer",
                "feature_id": "EH111111111",
                "associated_gene": gene_descriptors[6],
            },
            causative_event={
                "type": "CausativeEvent",
                "event_type": "read-through",
            },
            assay={
                "type": "Assay",
                "method_uri": "pmid:33576979",
                "assay_id": "obi:OBI_0003094",
                "assay_name": "fluorescence in-situ hybridization assay",
                "fusion_detection": "inferred",
            },
        )


def test_fusion_abstraction_validator(transcript_segments, linkers):
    """Test that instantiation of abstract fusion fails."""
    # can't create base fusion
    with pytest.raises(ValidationError) as exc_info:
        assert AbstractFusion(structural_elements=[transcript_segments[2], linkers[0]])
    check_validation_error(
        exc_info, "Value error, Cannot instantiate Fusion abstract class"
    )


def test_file_examples():
    """Test example JSON files."""
    # if this loads, then Pydantic validation was successful
    import fusor.examples as _  # noqa: F401


def test_model_examples():
    """Test example objects as provided in Pydantic config classes"""
    models = [
        FunctionalDomain,
        TranscriptSegmentElement,
        LinkerElement,
        TemplatedSequenceElement,
        GeneElement,
        UnknownGeneElement,
        MultiplePossibleGenesElement,
        RegulatoryElement,
        Assay,
        CausativeEvent,
        AssayedFusion,
        CategoricalFusion,
    ]
    for model in models:
        schema = model.model_config["json_schema_extra"]
        if "example" in schema:
            model(**schema["example"])
