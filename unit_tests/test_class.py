"""Unit testing for fusion class"""
import pytest
from fusion.model import GeneDescriptor, SequenceLocation, ChromosomeLocation
from fusion.model import GenomicRegion, TranscriptComponent, CriticalDomain
from fusion.model import Event, Linker, UnknownGene, RegulatoryElement, Fusion


def test_gene_descriptor():
    """Test GeneDescriptor object initializes correctly"""
    gen1 = GeneDescriptor(value_id='hgnc:1', label='G1')
    # value_id and label
    assert gen1.__dict__['value_id'] == 'hgnc:1'
    assert gen1.__dict__['label'] == 'G1'

    with pytest.raises(AssertionError):
        GeneDescriptor(value_id='hgnc2', label=5.3)


def test_sequence_location():
    """Test SequenceLocation object initializes correctly"""
    s1 = SequenceLocation(sequence_id='ncbi:NC_000001.11', start=1000,
                          end=5000)
    assert s1.__dict__['sequence_id'] == 'ncbi:NC_000001.11'
    assert s1.__dict__['start'] == 1000
    assert s1.__dict__['end'] == 5000

    with pytest.raises(AssertionError):
        SequenceLocation(sequence_id='ncbiNC_000001.11', start=1000.3,
                         end='5000')


def test_chromosome_location():
    """Test ChromosomeLocation object initializes correctly"""
    c1 = ChromosomeLocation(species='taxonomy:4232', chr='5', start='p13.2',
                            end='p13.5')
    assert c1.__dict__['species'] == 'taxonomy:4232'
    assert c1.__dict__['chr'] == '5'
    assert c1.__dict__['start'] == 'p13.2'
    assert c1.__dict__['end'] == 'p13.5'

    with pytest.raises(AssertionError):
        ChromosomeLocation(species='taxonomy4232', chr=5, start='s13.2',
                           end=14)


def test_event():
    """Test Event object initializes correctly"""
    event1 = Event(event_type='rearrangement')
    assert event1.__dict__['event_type'] == 'rearrangement'

    with pytest.raises(AssertionError):
        Event(event_type='combination')


def test_linker():
    """Test Linker object initializes correctly"""
    linker1 = Linker(linker_sequence='ATATGCA')
    assert linker1.__dict__['linker_sequence'] == 'ATATGCA'
    assert linker1.__dict__['linker_sequence'].lower() == 'atatgca'

    with pytest.raises(AssertionError):
        Linker(linker_sequence='ATDCHDGH')


def test_regulatory():
    """Test RegulatoryElement object initializes correctly"""
    r1 = RegulatoryElement(type='promoter', value_id='hgnc:435', label='G1')
    r2 = RegulatoryElement(type='enhancer', value_id='hgnc:435', label='G1')

    assert r1.__dict__['type'] == 'promoter'
    assert r1.__dict__['value_id'] == 'hgnc:435'
    assert r1.__dict__['label'] == 'G1'
    assert r2.__dict__['type'] == 'enhancer'
    assert r2.__dict__['value_id'] == 'hgnc:435'
    assert r2.__dict__['label'] == 'G1'

    with pytest.raises(AssertionError):
        RegulatoryElement(type='notpromoter', value_id='hgnc435', label=5)


def test_genomic_region_seq():
    """Test GenomicRegion object with SequenceLocation"""
    g1 = GenomicRegion(value=SequenceLocation(type='SequenceLocation',
                                              sequence_id='ncbi:NC_000001.11',
                                              start=154170399, end=1556346346))
    vals = g1.__dict__['value']
    assert vals.__dict__['type'] == 'SequenceLocation'
    assert vals.__dict__['sequence_id'] == 'ncbi:NC_000001.11'
    assert vals.__dict__['start'] == 154170399
    assert vals.__dict__['end'] == 1556346346

    with pytest.raises(AssertionError):
        GenomicRegion(value=SequenceLocation(type='GeneLocation',
                                             sequence_id='ncbiNC_000001.11',
                                             start='154170399',
                                             end='1556346346'))


def test_genomic_region_chr():
    """Test GenomicRegion object with ChromosomeLocation,"""
    g1 = GenomicRegion(value=ChromosomeLocation(type='ChromosomeLocation',
                                                species='taxonomy:9606',
                                                chr='12',
                                                start='p12.1', end='p12.2'))
    vals = g1.__dict__['value']
    assert vals.__dict__['type'] == 'ChromosomeLocation'
    assert vals.__dict__['species'] == 'taxonomy:9606'
    assert vals.__dict__['chr'] == '12'
    assert vals.__dict__['start'] == 'p12.1'
    assert vals.__dict__['end'] == 'p12.2'

    with pytest.raises(AssertionError):
        GenomicRegion(value=ChromosomeLocation(type='GeneLocation',
                                               species='taxonomy9606',
                                               chr=2, start='s12.1',
                                               end='s12.2'))


def test_transcript_component():
    """Test TranscriptComponent object initializes correctly"""
    gen1 = GeneDescriptor(value_id='hgnc:1', label='G1')
    g1 = GenomicRegion(value=ChromosomeLocation(type='ChromosomeLocation',
                                                species='taxonomy:9606',
                                                chr='12',
                                                start='p12.1', end='p12.2'))
    tr1 = TranscriptComponent(transcript='nm:152263.3', exon_start=1,
                              exon_start_offset=-9,
                              exon_end=8, exon_end_offset=7, gene=gen1,
                              component_genomic_region=g1)
    assert tr1.__dict__['transcript'] == 'nm:152263.3'
    assert tr1.__dict__['exon_start'] == 1
    assert tr1.__dict__['exon_start_offset'] == -9
    assert tr1.__dict__['exon_end'] == 8
    assert tr1.__dict__['exon_end_offset'] == 7
    gen = tr1.__dict__['gene']
    assert gen.__dict__['value_id'] == 'hgnc:1'
    assert gen.__dict__['label'] == 'G1'
    gr = tr1.__dict__['component_genomic_region'].__dict__['value']
    assert gr.__dict__['type'] == 'ChromosomeLocation'
    assert gr.__dict__['species'] == 'taxonomy:9606'
    assert gr.__dict__['chr'] == '12'
    assert gr.__dict__['start'] == 'p12.1'
    assert gr.__dict__['end'] == 'p12.2'

    with pytest.raises(AssertionError):
        TranscriptComponent(transcript='nm152263.3', exon_start='1',
                            exon_start_offset='-9',
                            exon_end='8', exon_end_offset='7', gene=4,
                            component_genomic_region=5)


def test_critical_domain():
    """Test CriticalDomain object initializes correctly"""
    gen1 = GeneDescriptor(value_id='hgnc:1', label='G1')
    domain = CriticalDomain(status='preserved',
                            name='tyrosine kinase catalytic domain',
                            id='interpro:IPR020635',
                            gene=gen1)

    assert domain.__dict__['status'] == 'preserved'
    assert domain.__dict__['name'] == 'tyrosine kinase catalytic domain'
    assert domain.__dict__['id'] == 'interpro:IPR020635'
    gen = domain.__dict__['gene']
    assert gen.__dict__['value_id'] == 'hgnc:1'
    assert gen.__dict__['label'] == 'G1'

    with pytest.raises(AssertionError):
        CriticalDomain(status='gained',
                       name='tyrosine kinase catalytic domain',
                       id='interproIPR020635',
                       gene=7)


def test_fusion():
    """Test Fusion Object initializes correctly"""
    gen1 = GeneDescriptor(value_id='hgnc:12012', label='TPM3')
    gen2 = GeneDescriptor(value_id='hgnc:8031', label='NTRK1')
    reg1 = RegulatoryElement(type='promoter', value_id='hgnc:100',
                             label='ASIC1')
    reg2 = RegulatoryElement(type='enhancer', value_id='hgnc:200', label='G1')

    domain = CriticalDomain(status='preserved',
                            name='tyrosine kinase catalytic domain',
                            id='interpro:IPR020635',
                            gene=gen2)

    ur1 = SequenceLocation(type='SequenceLocation',
                           sequence_id='ncbi:NC_000001.11',
                           start=1000, end=5000)
    c1 = UnknownGene(region=ur1)
    ur2 = ChromosomeLocation(type='ChromosomeLocation',
                             species='taxonomy:9606',
                             chr='12', start='p12.1', end='p12.16')
    c2 = UnknownGene(region=ur2)

    ex1 = GenomicRegion(value=SequenceLocation(type='SequenceLocation',
                                               sequence_id='ncbi:NC_000001.11',
                                               start=154170399,
                                               end=1556346346))
    ex2 = GenomicRegion(value=SequenceLocation(type='SequenceLocation',
                                               sequence_id='ncbi:NC_000001.13',
                                               start=154170401,
                                               end=1556346348))

    tr1 = TranscriptComponent(transcript='nm:152263.3', exon_start=1,
                              exon_start_offset=-9,
                              exon_end=8, exon_end_offset=7, gene=gen1,
                              component_genomic_region=ex1)
    tr2 = TranscriptComponent(transcript='nm:002529.3', exon_start=10,
                              exon_start_offset=-8,
                              exon_end=17, exon_end_offset=8, gene=gen2,
                              component_genomic_region=ex2)
    linker_seq = Linker(linker_sequence='ATTATTA')

    cause_eve = Event(event_type='rearrangement')

    # Validate correct input for Fusion object
    Fusion(r_frame_preserved=True, regulatory_elements=[reg1, reg2],
           protein_domains=[domain],
           transcript_components=[c1, c2, tr1, linker_seq, tr2],
           causative_event=cause_eve)

    with pytest.raises(AssertionError):
        Fusion(r_frame_preserved='T', regulatory_elements=[5, 5],
               protein_domains=[7, 7],
               transcript_components=[8], causative_event='Event')
