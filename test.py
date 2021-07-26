# pyGene Fusion local testing
import json
import re
from pydantic import BaseModel, validator, ValidationError, StrictStr, StrictInt, StrictBool
from typing import Optional, List, Union

class GeneDescriptor(BaseModel):
    """Define a Gene, denoted by GeneDescriptor"""
    type = 'GeneDescriptor'
    value_id: str
    label: StrictStr

    def __subscr__(self, attribute):
        return getattr(self, attribute)

    @validator('value_id')
    def is_valid(cls,v):
        assert v.count(':') == 1 and v.find(' ') == -1, 'id must be a CURIE'
        return v

class SequenceLocation(BaseModel):
    """Define SequenceLocation class"""
    type = 'SequenceLocation'
    seq_id: str
    start: StrictInt
    end: StrictInt

    @validator('seq_id')
    def is_curie(cls, v):
        assert v.count(':') == 1 and v.find(' ') == -1, 'chr must be a CURIE'
        return v

class ChromosomeLocation(BaseModel):
    """Define ChromosomeLocation class"""
    type = 'ChromosomeLocation'
    species: str
    chr: str
    start: str
    end: str

    @validator('species')
    def is_curie(cls, v):
        assert v.count(':') == 1 and v.find(' ') == -1, 'chr must be a CURIE'
        return v

    @validator('chr')
    def valid_chr(cls, v):
        assert v.isalnum(), 'chr must have only alphanumeric characters'
        return v

    @validator('start', 'end')
    def valid_loc(cls, v):
        assert bool(re.match('^cen|[pq](ter|([1-9][0-9]*(\.[1-9][0-9]*)?))$', v)) == True, \
            'start/end positions must match the regular expression ^cen|[pq](ter|([1-9][0-9]*(\.[1-9][0-9]*)?))$'
        return v

class GenomicRegion(BaseModel):
    """Define GenomicRegion class"""
    type = 'computed'
    description: Optional[str] = None
    value: Union[SequenceLocation, ChromosomeLocation]

class TranscriptComponent(BaseModel):
    """Define TranscriptComponent class, incorporating output from GeneDescriptor and GenomicRegion classes"""
    component_type = 'transcript_segment'
    transcript: str
    exon_start: StrictInt
    exon_start_offset : StrictInt = 0
    exon_end: StrictInt
    exon_end_offset : StrictInt = 0
    gene: GeneDescriptor
    component_genomic_region: GenomicRegion

    @validator('transcript')
    def correct_form(cls, v):
        assert v.count(':') == 1 and v.find(' ') == -1, 'id must be a CURIE'
        return v

class CriticalDomain(BaseModel):
    """Define CriticalDomain class"""
    status: StrictStr
    name: str
    id: str
    gene: GeneDescriptor

    @validator('status')
    def correct_status(cls, v):
        assert v.lower() == 'lost' or v == 'preserved', 'status must be either lost or preserved'
        return v

    @validator('id')
    def is_valid(cls, v):
        assert v.count(':') == 1 and v.find(' ') == -1, 'id must be a CURIE'
        return v

class Event(BaseModel):
    """Define Event class (causative event)"""
    event_type: str

    @validator('event_type')
    def event_validate(cls, v):
        assert v.lower() == 'rearrangement' or v == 'read-through' or v == 'trans-splicing', \
            'event entry must be one of rearrangement, read-through, or trans-splicing'
        return v

class Linker(BaseModel):
    """Define Linker class (linker sequence)"""
    linker_sequence: str
    component_type = 'linker_sequence'

    @validator('linker_sequence')
    def valid_input(cls, v):
        assert set(v.upper()) <= set('ACGT'), 'Linker sequence must only contain A,C,G,T'
        return v

class UnknownGene(BaseModel):
    """Define UnknownGene class"""
    component_type = 'unknown_gene'
    region: Optional[Union[SequenceLocation, ChromosomeLocation]]

class RegulatoryElement(BaseModel):
    """Define RegulatoryElement class"""
    type: str
    value_id: str
    label: StrictStr

    @validator('type')
    def valid_reg_type(cls, v):
        assert v.lower() == 'promoter' or v == 'enhancer', 'type must be either promoter or enhancer'
        return v

    @validator('value_id')
    def valid_id(cls, v):
        assert v.count(':') == 1 and v.find(' ') == -1, 'id must be a CURIE'
        return v

class Fusion(BaseModel):
    """Define Fusion class"""
    r_frame_preserved: StrictBool
    regulatory_elements: List[RegulatoryElement]
    protein_domains: List[CriticalDomain]
    transcript_components: List[Union[TranscriptComponent, Linker, UnknownGene]]
    causative_event: Event

    def make_json(self):
        return(json.dumps(self, default=lambda o: o.__dict__, indent=4))

def main():
    try:
        gen1 = GeneDescriptor(value_id = 'hgnc:12012', label = 'TPM3')
        gen2 = GeneDescriptor(value_id = 'hgnc:8031', label = 'NTRK1')
        ex1 = GenomicRegion(value = SequenceLocation(type = 'SequenceLocation', seq_id = 'ncbi:NC_000001.11',  start =  154170399, end = 1556346346))
        ex2 = GenomicRegion(value = SequenceLocation(type = 'SequenceLocation', seq_id = 'ncbi:NC_000001.13',  start =  154170401, end = 1556346348))
        linker_seq = Linker(linker_sequence = 'ATTATTA')
        cause_eve = Event(event_type = 'rearrangement')
        domain = CriticalDomain(status = 'preserved', name = 'tyrosine kinase catalytic domain', id = 'interpro:IPR020635',
                         gene = gen2)
        reg1 = RegulatoryElement(type = 'promoter', value_id = 'hgnc:100', label = 'ASIC1')
        reg2 = RegulatoryElement(type = 'enhancer', value_id = 'hgnc:200', label = 'G1')
        ur1 = SequenceLocation(type = 'SequenceLocation', seq_id = 'ncbi:NC_000001.11', start = 1000, end = 5000)
        c1 = UnknownGene(region = ur1)
        c2 = UnknownGene()
        ur2 = ChromosomeLocation(type = 'ChromosomeLocation', species = 'taxonomy:9606' ,chr = '12', start = 'p12.1', end = 'p12.16')
        c3 = UnknownGene(region = ur2)
        tr1 = TranscriptComponent(transcript = 'nm:152263.3', exon_start = 1, exon_start_offset = -9,
                                  exon_end = 8, exon_end_offset = 7, gene = gen1, component_genomic_region = ex1)
        tr2 = TranscriptComponent(transcript = 'nm:002529.3', exon_start = 10, exon_start_offset = -8,
                                  exon_end = 17, exon_end_offset = 8, gene = gen2, component_genomic_region = ex2)
        molecular_fusion = Fusion(r_frame_preserved = True, regulatory_elements = [reg1, reg2], protein_domains = [domain],
                                  transcript_components = [c1, c2, c3, tr1, linker_seq, tr2], causative_event = cause_eve)
        print(molecular_fusion.json())
    except ValidationError as f:
        print(f)

    filename = 'example_gene_fusion_class.json'
    with open(filename, 'w') as files:
        files.write(molecular_fusion.make_json())
    files.close()


if __name__ == '__main__':
    main()

