"""Module for testing the fusion model."""
import pydantic
import pytest
from fusor.model import GeneDescriptor, SequenceLocation, ChromosomeLocation,\
    GenomicRegion, TranscriptComponent, CriticalDomain, Event, Linker, \
    UnknownGene, RegulatoryElement, Fusion, Extension, GeneValueObject,\
    SimpleInterval, CytobandInterval


def test_simpleinterval():
    """Test SimpleInterval object"""
    s1 = SimpleInterval(start=2, end=5)
    assert s1.__dict__['start'] == 2
    assert s1.__dict__['end'] == 5

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        SimpleInterval(start='s', end='q')


def test_cytobandinterval():
    """Test CytobandInterval object"""
    c1 = CytobandInterval(start='q17.1', end='q17.2')
    assert c1.__dict__['start'] == 'q17.1'
    assert c1.__dict__['end'] == 'q17.2'

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        CytobandInterval(start=2, end=5)


def test_genevalueobject():
    """Test GeneValueObject"""
    gen1 = GeneValueObject(id='hgnc:1')
    assert gen1.__dict__['id'] == 'hgnc:1'

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        GeneValueObject(id='hgnc1')


def test_extension():
    """Test Extension"""
    ex1 = Extension(name='strand', value='-')
    assert ex1.__dict__['name'] == 'strand'
    assert ex1.__dict__['value'] == '-'

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        Extension(name=3, value='-')


def test_gene_descriptor():
    """Test GeneDescriptor object initializes correctly"""
    gen1 = GeneDescriptor(id='test:g1', value=GeneValueObject(id='hgnc:1497',
                          type='Gene'), label='G1', xrefs=['ncbi:g1'],
                          alternate_labels=['gen1'])

    assert gen1.__dict__['id'] == 'test:g1'
    vals = gen1.__dict__['value']
    assert vals.__dict__['id'] == 'hgnc:1497'
    assert vals.__dict__['type'] == 'Gene'
    assert gen1.__dict__['label'] == 'G1'

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        GeneDescriptor(id=1, value=GeneValueObject(id='hgnc1497',
                       type='Gene'), label='G1', xrefs=['ncbi:g1'],
                       alternate_labels=['gen1'])


def test_sequence_location():
    """Test SequenceLocation object initializes correctly"""
    s1 = SequenceLocation(sequence_id='ncbi:NC_000001.11',
                          interval=SimpleInterval(start=22525, end=252525),
                          type='SequenceLocation')
    assert s1.__dict__['sequence_id'] == 'ncbi:NC_000001.11'
    si = s1.__dict__['interval']
    assert si.__dict__['start'] == 22525
    assert si.__dict__['end'] == 252525

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        SequenceLocation(sequence_id='ncbiNC_000001.11',
                         interval=SimpleInterval(start='s', end='p'))


def test_chromosome_location():
    """Test ChromosomeLocation object initializes correctly"""
    c1 = ChromosomeLocation(species_id='taxonomy:4232', chr='5',
                            interval=CytobandInterval(start='p5.2',
                                                      end='p5.3'),
                            type='ChromosomeLocation')
    assert c1.__dict__['species_id'] == 'taxonomy:4232'
    assert c1.__dict__['chr'] == '5'
    ci = c1.__dict__['interval']
    assert ci.__dict__['start'] == 'p5.2'
    assert ci.__dict__['end'] == 'p5.3'

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        ChromosomeLocation(species_ud='taxonomy:4232', chr='5',
                           interval=CytobandInterval(start=17,
                                                     end=18))


def test_event():
    """Test Event object initializes correctly"""
    event1 = Event(event_type='rearrangement')
    assert event1.__dict__['event_type'] == 'rearrangement'

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        Event(event_type='combination')


def test_linker():
    """Test Linker object initializes correctly"""
    linker1 = Linker(linker_sequence='ATATGCA')
    assert linker1.__dict__['linker_sequence'] == 'ATATGCA'
    assert linker1.__dict__['linker_sequence'].lower() == 'atatgca'

    with pytest.raises(pydantic.error_wrappers.ValidationError):
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

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        RegulatoryElement(type='notpromoter', value_id='hgnc435', label=5)


def test_genomic_region_seq():
    """Test GenomicRegion object with SequenceLocation"""
    si = SimpleInterval(start=15455, end=15566)
    g1 = GenomicRegion(value=SequenceLocation(sequence_id='ncbi:NC_000001.11',
                                              interval=si,
                                              type='SequenceLocation'))

    vals = g1.__dict__['value']
    assert vals.__dict__['sequence_id'] == 'ncbi:NC_000001.11'
    si = vals.__dict__['interval']
    assert si.__dict__['start'] == 15455
    assert si.__dict__['end'] == 15566

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        GenomicRegion(value=SequenceLocation(sequence_id='ncbiNC_000001.11',
                                             interval=si,
                                             type='SeqLocation'))


def test_genomic_region_chr():
    """Test GenomicRegion object with ChromosomeLocation,"""
    ci = CytobandInterval(start='p12.1', end='p12.2')
    g1 = GenomicRegion(value=ChromosomeLocation(species_id='taxonomy:9606',
                                                chr='12',
                                                interval=ci,
                                                type='ChromosomeLocation'))

    vals = g1.__dict__['value']
    assert vals.__dict__['species_id'] == 'taxonomy:9606'
    assert vals.__dict__['chr'] == '12'
    ci = vals.__dict__['interval']
    assert ci.__dict__['start'] == 'p12.1'
    assert ci.__dict__['end'] == 'p12.2'

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        GenomicRegion(value=ChromosomeLocation(species='taxonomy9606',
                                               chr=12,
                                               interval=ci,
                                               type='SequenceLocation'))


def test_transcript_component():
    """Test TranscriptComponent object initializes correctly"""
    gen1 = GeneDescriptor(id='test:1', value=GeneValueObject(id='hgnc:1'),
                          label='G1')
    ci = CytobandInterval(start='p12.1', end='p12.2')
    g1 = GenomicRegion(value=ChromosomeLocation(species_id='taxonomy:9606',
                                                chr='12',
                                                interval=ci,
                                                type='ChromosomeLocation'))
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
    assert gen.__dict__['id'] == 'test:1'
    assert gen.__dict__['label'] == 'G1'
    gv = gen.__dict__['value']
    assert gv.__dict__['id'] == 'hgnc:1'
    gr = tr1.__dict__['component_genomic_region'].__dict__['value']
    assert gr.__dict__['species_id'] == 'taxonomy:9606'
    assert gr.__dict__['type'] == 'ChromosomeLocation'
    assert gr.__dict__['chr'] == '12'
    grc = gr.__dict__['interval']
    assert grc.__dict__['start'] == 'p12.1'
    assert grc.__dict__['end'] == 'p12.2'

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        TranscriptComponent(transcript='nm152263.3', exon_start='1',
                            exon_start_offset='-9',
                            exon_end='8', exon_end_offset='7', gene=4,
                            component_genomic_region=5)


def test_critical_domain():
    """Test CriticalDomain object initializes correctly"""
    gen1 = GeneDescriptor(id='test:1', value=GeneValueObject(id='hgnc:1'),
                          label='G1')
    domain = CriticalDomain(status='preserved',
                            name='tyrosine kinase catalytic domain',
                            id='interpro:IPR020635',
                            gene=gen1)

    assert domain.__dict__['status'] == 'preserved'
    assert domain.__dict__['name'] == 'tyrosine kinase catalytic domain'
    assert domain.__dict__['id'] == 'interpro:IPR020635'
    gen = domain.__dict__['gene']
    assert gen.__dict__['id'] == 'test:1'
    assert gen.__dict__['label'] == 'G1'
    gv = gen.__dict__['value']
    assert gv.__dict__['id'] == 'hgnc:1'

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        CriticalDomain(status='gained',
                       name='tyrosine kinase catalytic domain',
                       id='interproIPR020635',
                       gene=7)


def test_fusion():
    """Test Fusion Object initializes correctly"""
    gen1 = GeneDescriptor(id='test:1', value=GeneValueObject(id='hgnc:1'),
                          label='G1')
    gen2 = GeneDescriptor(id='test:2', value=GeneValueObject(id='hgnc:2'),
                          label='G2', xrefs=['ncbigene:6', 'ensembl:7'])
    reg1 = RegulatoryElement(type='promoter', value_id='hgnc:100',
                             label='ASIC1')
    reg2 = RegulatoryElement(type='enhancer', value_id='hgnc:200', label='G1')

    domain = CriticalDomain(status='preserved',
                            name='tyrosine kinase catalytic domain',
                            id='interpro:IPR020635',
                            gene=gen2)

    ur1 = SequenceLocation(type='SequenceLocation',
                           sequence_id='ncbi:NC_000001.11',
                           interval=SimpleInterval(start=1000, end=5000))
    c1 = UnknownGene(region=ur1)
    cex = CytobandInterval(start='p12.1', end='p12.16')
    ur2 = ChromosomeLocation(type='ChromosomeLocation',
                             species_id='taxonomy:9606',
                             chr='12',
                             interval=cex)
    c2 = UnknownGene(region=ur2)
    si = SimpleInterval(start=15455, end=15566)
    ex1 = GenomicRegion(value=SequenceLocation(type='SequenceLocation',
                                               sequence_id='ncbi:NC_000001.11',
                                               interval=si))
    si2 = SimpleInterval(start=256235, end=296969)
    ex2 = GenomicRegion(value=SequenceLocation(type='SequenceLocation',
                                               sequence_id='ncbi:NC_000001.13',
                                               interval=si2))

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

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        Fusion(r_frame_preserved='T', regulatory_elements=[5, 5],
               protein_domains=[7, 7],
               transcript_components=[8], causative_event='Event')
