# pyGene Fusion local testing
import json
from pydantic import BaseModel, validator, ValidationError, StrictStr, StrictInt, StrictBool
from typing import Optional, List, Union

class Gene(BaseModel):
    type = 'GeneDescriptor'
    value_id: str
    label: StrictStr

    @validator('value_id')
    def is_valid(cls,v):
        assert v.find(':') != -1 and v.find(' ') == -1, 'id must be a CURIE'
        return v

class SequenceRegion(BaseModel):
    type = 'SequenceLocation'
    seq_id: str
    start: StrictInt
    end: StrictInt

    @validator('seq_id')
    def is_curie(cls, v):
        assert v.find(':') != -1 and v.find(' ') == -1, 'chr must be a CURIE'
        return v

class ChromosomeRegion(BaseModel):
    type = 'ChromosomeLocation'
    chr: str
    start: str
    end: str

    @validator('chr')
    def valid_chr(cls, v):
        assert v.isalnum(), 'chr must have only alphanumeric characters'
        return v

    @validator('start', 'end')
    def valid_loc(cls, v):
        assert v[0] == 'p' or v[0] == 'q', 'start/end must start with p or q'
        return v

class GenomicRegion(BaseModel):
    type = 'computed'
    description: Optional[str] = None
    value: Union[SequenceRegion, ChromosomeRegion]

class TranscriptRegion(BaseModel):
    component_type = 'transcript_segment'
    transcript: str
    exon_start: StrictInt
    exon_start_offset : StrictInt
    exon_end: StrictInt
    exon_end_offset : StrictInt
    gene: Gene
    component_genomic_region: GenomicRegion

    @validator('transcript')
    def correct_form(cls, v):
        assert v.find(':') != -1 and v.find(' ') == -1, 'id must be a CURIE'
        return v

class CriticalDomain(BaseModel):
    status: StrictStr
    name: str
    id: str
    gene: Gene

    @validator('status')
    def correct_status(cls, v):
        assert v == 'lost' or v == 'preserved', 'status must be either lost or preserved'
        return v

    @validator('id')
    def is_valid(cls, v):
        assert v.find(':') != -1 and v.find(' ') == -1, 'id must be a CURIE'
        return v

class Event(BaseModel):
    event_type: str

    @validator('event_type')
    def event_validate(cls, v):
        assert v == 'rearrangement' or v == 'read-through' or v == 'trans-splicing', \
            'event entry must be one of rearrangement, read-through, or trans-splicing'
        return v

class Linker(BaseModel):
    linker_sequence: str

    @validator('linker_sequence')
    def valid_input(cls, v):
        assert set(v) <= set('ACGT'), 'Linker sequence must only contain A,C,G,T'
        return v


class UnknownComponent(BaseModel):
    component_type = 'unknown_gene'
    region: Optional[Union[SequenceRegion, ChromosomeRegion]]

class RegulatoryElement(BaseModel):
    type: str
    value_id: str
    label: StrictStr

    @validator('type')
    def valid_reg_type(cls, v):
        assert v == 'promoter' or v == 'enhancer', 'type must be either promoter or enhancer'
        return v

    @validator('value_id')
    def valid_id(cls, v):
        assert v.find(':') != -1 and v.find(' ') == -1, 'id must be a CURIE'
        return v

class Fusion(BaseModel):
    r_frame_preserved: StrictBool
    regulatory_elements: List[RegulatoryElement]
    protein_domains: List[CriticalDomain]
    transcript_components = []
    causative_event: Event

    def make_json(self):
        return(json.dumps(self, default=lambda o: o.__dict__, indent=4))

def main():

    try:
        gen1 = Gene(value_id = 'hgnc:12012', label = 'TPM3')
        gen2 = Gene(value_id = 'hgnc:8031', label = 'NTRK1')
        ex1 = GenomicRegion(value = SequenceRegion(type = 'SequenceLocation', seq_id = 'ncbi:NC_000001.11',  start =  154170399, end = 1556346346))
        ex2 = GenomicRegion(value = SequenceRegion(type = 'SequenceLocation', seq_id = 'ncbi:NC_000001.13',  start =  154170401, end = 1556346348))
        linker_seq = Linker(linker_sequence = 'ATTATTA')
        cause_eve = Event(event_type = 'rearrangement')
        domain = CriticalDomain(status = 'preserved', name = 'tyrosine kinase catalytic domain', id = 'interpro:IPR020635',
                         gene = gen2)
        reg1 = RegulatoryElement(type = 'promoter', value_id = 'hgnc:100', label = 'ASIC1')
        reg2 = RegulatoryElement(type = 'enhancer', value_id = 'hgnc:200', label = 'G1')
        ur1 = SequenceRegion(type = 'SequenceLocation', seq_id = 'ncbi:NC_000001.11', start = 1000, end = 5000)
        c1 = UnknownComponent(region = ur1)
        c2 = UnknownComponent()
        ur2 = ChromosomeRegion(type = 'ChromosomeLocation', chr = '12', start = 'p12.1', end = 'p12.16')
        c3 = UnknownComponent(region = ur2)
        tr1 = TranscriptRegion(transcript = 'NM:152263.3', exon_start = 1, exon_start_offset = -9,
                                  exon_end = 8, exon_end_offset = 7, gene = gen1, component_genomic_region = ex1)
        tr2 = TranscriptRegion(transcript = 'NM:002529.3', exon_start = 10, exon_start_offset = -8,
                                  exon_end = 17, exon_end_offset = 8, gene = gen2, component_genomic_region = ex2)
        molecular_fusion = Fusion(r_frame_preserved = True, regulatory_elements = [reg1, reg2], protein_domains = [domain],
                                  transcript_components = [c1, c2, c3, tr1, linker_seq, tr2], causative_event = cause_eve,)
        print(molecular_fusion.json())
    except ValidationError as f:
        print(f)

    filename = 'example_gene_fusion_class.json'
    with open(filename, 'w') as files:
        files.write(molecular_fusion.make_json())
    files.close()


if __name__ == '__main__':
    main()

