"""Module for testing the fusion model."""
from pydantic import ValidationError
import pytest
from fusor.model import LocationDescriptor, TranscriptSegmentComponent, \
    GenomicRegionComponent, UnknownGeneComponent, GeneComponent, \
    LinkerComponent, SequenceDescriptor, CriticalDomain, Event, \
    RegulatoryElement, Fusion


@pytest.fixture(scope='module')
def gene_descriptors():
    """Provide possible gene_descriptor input."""
    return [
        {
            'id': 'gene:G1',
            'value': {'id': 'hgnc:9339'},
            'label': 'G1'
        },
        {
            'id': 'gene:ABL',
            'value': {'id': 'hgnc:76'},
            'label': 'ABL'
        },
        {
            'id': 'gene:BCR1',
            'value': {'id': 'hgnc:1014'},
            'label': 'BCR1'
        },
        {
            'id': 'gene:NTRK1',
            'value_id': 'hgnc:8031',
            'label': 'NTRK1'
        }
    ]


@pytest.fixture(scope='module')
def critical_domains(gene_descriptors):
    """Provide possible critical_domains input."""
    return [
        {
            'status': 'preserved',
            'name': 'phlebovirus glycoprotein g1',
            'id': 'interpro:IPR010826',
            'gene': gene_descriptors[0]
        },
        {
            'status': 'lost',
            'name': 'tyrosine kinase catalytic domain',
            'id': 'interpro:IPR020635',
            'gene': gene_descriptors[3]
        }
    ]


@pytest.fixture(scope='module')
def location_descriptors():
    """Provide possible genomic_region input."""
    return [
        {
            'id': 'NC_000001.11:15455-15566',
            'type': 'LocationDescriptor',
            'value': {
                'sequence_id': 'ncbi:NC_000001.11',
                'interval': {'start': 15455, 'end': 15566},
                'type': 'SequenceLocation'
            },
            'label': 'NC_000001.11:15455-15566',
        },
        {
            'id': 'chr12:p12.1-p12.2',
            'type': 'LocationDescriptor',
            'value': {
                'species_id': 'taxonomy:9606',
                'chr': '12',
                'interval': {'start': 'p12.1', 'end': 'p12.2'}
            },
            'label': 'chr12:p12.1-p12.2',
        },
    ]


@pytest.fixture(scope='module')
def transcript_segments(location_descriptors):
    """Provide possible transcript_segment input."""
    return [
        {
            'transcript': 'refseq:NM_152263.3',
            'exon_start': 1,
            'exon_start_offset': -9,
            'exon_end': 8,
            'exon_end_offset': 7,
            'gene': {'id': 'test:1', 'value': {'id': 'hgnc:1'}, 'label': 'G1'},
            'component_genomic_region': location_descriptors[1]
        },
        {
            'component_type': 'transcript_segment',
            'transcript': 'refseq:NM_034348.3',
            'exon_start': 1,
            'exon_end': 8,
            'gene': {
                'id': 'gene:NTRK1',
                'value_id': 'hgnc:8031',
                'label': 'NTRK1',
            },
            'component_genomic_region': location_descriptors[0]
        },
        {
            'component_type': 'transcript_segment',
            'transcript': 'refseq:NM_938439.4',
            'exon_start': 7,
            'exon_end': 14,
            'exon_end_offset': -5,
            'gene': {
                'id': 'gene:ALK',
                'value_id': 'hgnc:1837',
                'label': 'ALK',
            },
            'component_genomic_region': location_descriptors[0]
        }
    ]


@pytest.fixture(scope='module')
def gene_components():
    """Provide possible gene component input data."""
    return [
        {
            'component_type': 'gene',
            'gene': {
                'id': 'gene:ABL1',
                'value_id': 'hgnc:2000',
            }
        }
    ]


@pytest.fixture(scope='module')
def genomic_region_components(location_descriptors):
    """Provide possible genomic_region component input data."""
    return [
        {
            'component_type': 'genomic_region',
            'strand': '+',
            'region': location_descriptors[1]
        },
        {
            'component_type': 'genomic_region',
            'strand': '-',
            'region': location_descriptors[0]
        }
    ]


@pytest.fixture(scope='module')
def sequence_descriptors():
    """Provide possible SequenceDescriptor input data"""
    return [
        {
            'id': 'sequence:ACGT',
            'type': 'SequenceDescriptor',
            'value': {
                'sequence': 'ACGT',
                'type': 'SequenceState',
            },
            'residue_type': 'SO:0000348'
        },
        {
            'id': 'sequence:ATGTA',
            'type': 'SequenceDescriptor',
            'value': {
                'sequence': 'ATGTA',
                'type': 'SequenceState',
            },
            'residue_type': 'SO:0000348'
        }
    ]


@pytest.fixture(scope='module')
def linkers(sequence_descriptors):
    """Provide possible linker component input data."""
    return [
        {
            'component_type': 'linker_sequence',
            'linker_sequence': sequence_descriptors[0]
        },
        {
            'component_type': 'linker_sequence',
            'linker_sequence': sequence_descriptors[1]
        }
    ]


@pytest.fixture(scope='module')
def regulatory_elements(gene_descriptors):
    """Provide possible regulatory_element input data."""
    return [
        {
            'type': 'promoter',
            'gene': gene_descriptors[0]
        }
    ]


def test_critical_domain(critical_domains, gene_descriptors):
    """Test CriticalDomain object initializes correctly"""
    test_domain = CriticalDomain(**critical_domains[0])
    assert test_domain.status == 'preserved'
    assert test_domain.name == 'phlebovirus glycoprotein g1'
    assert test_domain.id == 'interpro:IPR010826'
    assert test_domain.gene.id == 'gene:G1'
    assert test_domain.gene.label == 'G1'
    assert test_domain.gene.value.id == 'hgnc:9339'

    # test status string
    with pytest.raises(ValidationError):
        CriticalDomain(**{
            'status': 'gained',
            'name': 'tyrosine kinase catalytic domain',
            'id': 'interpro:IPR020635',
            'gene': gene_descriptors[0]
        })

    # test domain ID CURIE requirement
    with pytest.raises(ValidationError):
        CriticalDomain(**{
            'status': 'lost',
            'name': 'tyrosine kinase catalytic domain',
            'id': 'interpro_IPR020635',
            'gene': gene_descriptors[0]
        })


def test_location_descriptors(location_descriptors):
    """Test GenomicRegion class."""
    # check with SequenceLocation
    test_loc_descr = LocationDescriptor(**location_descriptors[0])
    assert test_loc_descr.id == 'NC_000001.11:15455-15566'
    assert test_loc_descr.type == 'LocationDescriptor'
    assert test_loc_descr.value.sequence_id == 'ncbi:NC_000001.11'
    assert test_loc_descr.value.interval.start == 15455
    assert test_loc_descr.value.interval.end == 15566
    assert test_loc_descr.value.type == 'SequenceLocation'
    assert test_loc_descr.label == 'NC_000001.11:15455-15566'

    # ID required
    with pytest.raises(ValidationError):
        LocationDescriptor(**{
            'value': {
                'sequence_id': 'ncbi:NC_000001.11',
                'interval': {
                    'type': 'SimpleInterval',
                    'end': 23490823409,
                    'start': 23940834,
                },
                'type': 'SequenceLocation'
            },
            'type': 'LocationDescriptor'
        })

    # check with ChromosomeLocation
    test_loc_descr = LocationDescriptor(**location_descriptors[1])
    assert test_loc_descr.id == 'chr12:p12.1-p12.2'
    assert test_loc_descr.type == 'LocationDescriptor'
    assert test_loc_descr.value.species_id == 'taxonomy:9606'
    assert test_loc_descr.value.chr == '12'
    assert test_loc_descr.value.interval.start == 'p12.1'
    assert test_loc_descr.value.interval.end == 'p12.2'
    assert test_loc_descr.label == 'chr12:p12.1-p12.2'

    # ID must be CURIE
    with pytest.raises(ValidationError):
        LocationDescriptor(**{
            'id': 'chr12-p12.1-p12.2',
            'value': {
                'species_id': 'taxonomy:9606',
                'chr': '12',
                'interval': {
                    'start': 'p12.1', 'end': 'p12.2',
                    'type': 'CytobandInterval'
                },
                'type': 'ChromosomeLocation'
            },
        })


def test_transcript_segment_component(transcript_segments):
    """Test TranscriptSegmentComponent object initializes correctly"""
    test_component = TranscriptSegmentComponent(**transcript_segments[0])
    assert test_component.transcript == 'refseq:NM_152263.3'
    assert test_component.exon_start == 1
    assert test_component.exon_start_offset == -9
    assert test_component.exon_end == 8
    assert test_component.exon_end_offset == 7
    assert test_component.gene.id == 'test:1'
    assert test_component.gene.label == 'G1'
    assert test_component.gene.value.id == 'hgnc:1'
    test_region = test_component.component_genomic_region
    assert test_region.value.species_id == 'taxonomy:9606'
    assert test_region.value.type == 'ChromosomeLocation'
    assert test_region.value.chr == '12'
    assert test_region.value.interval.start == 'p12.1'
    assert test_region.value.interval.end == 'p12.2'

    # check CURIE requirement
    with pytest.raises(ValidationError):
        TranscriptSegmentComponent(**{
            'transcript': 'NM_152263.3',
            'exon_start': '1',
            'exon_start_offset': '-9',
            'exon_end': '8',
            'exon_end_offset': '7',
            'gene': {'id': 'test:1', 'value': {'id': 'hgnc:1'}, 'label': 'G1'},
            'component_genomic_region': {
                'value': {
                    'species_id': 'taxonomy:9606', 'chr': '12',
                    'interval': {'start': 'p12.1', 'end': 'p12.2'},
                }
            }
        })


def test_sequence_descriptor(sequence_descriptors):
    """Test that SequenceDescriptor objects initialize correctly"""
    test_descriptor = SequenceDescriptor(**sequence_descriptors[0])
    assert test_descriptor.id == 'sequence:ACGT'
    assert test_descriptor.type == 'SequenceDescriptor'
    assert test_descriptor.value.sequence == 'ACGT'
    assert test_descriptor.value.type == 'SequenceState'
    assert test_descriptor.residue_type == 'SO:0000348'


def test_linker_component(linkers):
    """Test Linker object initializes correctly"""
    test_linker = LinkerComponent(**linkers[0])
    assert test_linker.component_type == 'linker_sequence'
    assert test_linker.linker_sequence.id == 'sequence:ACGT'
    assert test_linker.linker_sequence.type == 'SequenceDescriptor'
    assert test_linker.linker_sequence.value.sequence == 'ACGT'
    assert test_linker.linker_sequence.value.type == 'SequenceState'
    assert test_linker.linker_sequence.residue_type == 'SO:0000348'

    # check base validation
    with pytest.raises(ValidationError):
        LinkerComponent(**{
            'linker_sequence':
                {
                    'id': 'sequence:ABGT',
                    'value': {
                        'sequence': 'ABGT'
                    }
                }
        })


def test_genomic_region_component(genomic_region_components):
    """Test that GenomicRegionComponent initializes correctly."""
    test_component = GenomicRegionComponent(**genomic_region_components[0])
    assert test_component.component_type == 'genomic_region'
    assert test_component.strand.value == '+'
    assert test_component.region.id == 'chr12:p12.1-p12.2'
    assert test_component.region.type == 'LocationDescriptor'
    assert test_component.region.value.species_id == 'taxonomy:9606'
    assert test_component.region.value.chr == '12'
    assert test_component.region.value.interval.start == 'p12.1'
    assert test_component.region.value.interval.end == 'p12.2'
    assert test_component.region.label == 'chr12:p12.1-p12.2'

    with pytest.raises(ValidationError):
        GenomicRegionComponent(**{
            'region': {
                'interval': {
                    'start': 39408,
                    'stop': 39414
                }
            },
            'sequence_id': 'ga4gh:SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl'
        })


def test_gene_component(gene_descriptors):
    """Test that Gene component initializes correctly."""
    test_component = GeneComponent(**{'gene': gene_descriptors[0]})
    assert test_component.component_type == 'gene'
    assert test_component.gene.id == 'gene:G1'
    assert test_component.gene.label == 'G1'
    assert test_component.gene.value.id == 'hgnc:9339'

    # test CURIE requirement
    with pytest.raises(ValidationError):
        GenomicRegionComponent(**{
            'gene': {
                'id': 'G1',
                'value': {'id': 'hgnc:9339'},
                'label': 'G1'
            }
        })


def test_unknown_gene_component():
    """Test that unknown_gene component initializes correctly."""
    test_component = UnknownGeneComponent()
    assert test_component.component_type == 'unknown_gene'


def test_event():
    """Test Event object initializes correctly"""
    rearrangement = 'rearrangement'
    test_event = Event(rearrangement)
    assert test_event.value == rearrangement

    with pytest.raises(ValueError):
        Event('combination')


def test_regulatory(regulatory_elements, gene_descriptors):
    """Test RegulatoryElement object initializes correctly"""
    test_reg_elmt = RegulatoryElement(**regulatory_elements[0])
    assert test_reg_elmt.type.value == 'promoter'
    assert test_reg_elmt.gene.id == 'gene:G1'
    assert test_reg_elmt.gene.value.id == 'hgnc:9339'
    assert test_reg_elmt.gene.label == 'G1'

    # check type constraint
    with pytest.raises(ValidationError):
        RegulatoryElement(**{
            'type': 'notpromoter',
            'gene': gene_descriptors[0]
        })


def test_fusion(critical_domains, transcript_segments,
                genomic_region_components, linkers, gene_components,
                regulatory_elements):
    """Test that Fusion object initializes correctly"""
    unknown_component = {
        'component_type': 'unknown_gene',
    }

    # test valid object
    fusion = Fusion(**{
        'r_frame_preserved': True,
        'protein_domains': [critical_domains[0]],
        'transcript_components': [
            transcript_segments[1], transcript_segments[2]
        ],
        'causative_event': 'rearrangement',
        'regulatory_elements': [regulatory_elements[0]]
    })

    assert fusion.transcript_components[0].transcript == 'refseq:NM_034348.3'

    # test that non-component properties are optional
    assert Fusion(**{
        'transcript_components': [
            transcript_segments[1], transcript_segments[2]
        ]
    })

    # test variety of component types
    assert Fusion(**{
        'transcript_components': [
            unknown_component,
            gene_components[0],
            transcript_segments[2],
            linkers[1],
            genomic_region_components[1]
        ]
    })

    # components are mandatory
    with pytest.raises(ValidationError):
        assert Fusion(**{
            'r_frame_preserved': True,
            'protein_domains': [critical_domains[1]],
            'causative_event': 'rearrangement',
            'regulatory_elements': [regulatory_elements[0]]
        })

    # must have >= 2 components
    with pytest.raises(ValidationError):
        assert Fusion(**{
            'transcript_components': [unknown_component]
        })
