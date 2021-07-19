# pyGene Fusion local testing
import json
from pydantic import BaseModel, validator, ValidationError, StrictStr, StrictInt, StrictBool

class Gene(BaseModel):
    id: str
    symbol : StrictStr

    @validator('id')
    def is_valid(cls,v):
        assert v.find(':') != -1 and v.find(' ') == -1, 'id must be a CURIE'
        return v

class Exon_loc(BaseModel):
    chr : str
    position: StrictInt

    @validator('chr')
    def is_valid_chr(cls, v):
        assert v[0:3] == 'NC_', 'chr must start with NC_'
        return v

class Transcript_region(BaseModel):
    transcript : str
    gene : Gene
    exon_start : StrictInt
    exon_start_genomic : Exon_loc
    exon_end : StrictInt
    exon_end_genomic: Exon_loc

    @validator('transcript')
    def correct_form(cls, v):
        assert v[0:3] == 'NM_', 'transcript must start with NM_'
        return v

class Domains(BaseModel):
    status : StrictStr
    name: str
    id : str
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
    linker_sequence : str

    @validator('linker_sequence')
    def valid_input(cls, v):
        assert set(v) <= set('ACGT'), 'Linker sequence must only contain A,C,G,T'
        return v

class Fusion(BaseModel):
    r_frame_preserved: StrictBool
    domains : Domains
    components = []
    causative_event: Event

    def make_json(self):
        return(json.dumps(self, default=lambda o: o.__dict__, indent=4))

def main():

    try:
        gen1 = Gene(id = 'hgnc:12012', symbol = 'TPM3')
        gen2 = Gene(id = 'hgnc:8031', symbol = 'NTRK1')
        ex1 = Exon_loc(chr = 'NC_000001.11', position =  154170399)
        ex2 = Exon_loc(chr = 'NC_000001.11', position = 158843425)
        ex3 = Exon_loc(chr = 'NC_0000011.1', position = 154545)
        ex4 = Exon_loc(chr = 'NC_45634534', position = 68686868)
        linker_seq = Linker(linker_sequence = 'ATTATTA')
        cause_eve = Event(event_type = 'rearrangement')
        domain = Domains(status = 'preserved', name = 'tyrosine kinase catalytic domain', id = 'interpro:IPR020635',
                         gene = gen2)
        tr1 = Transcript_region(transcript = 'NM_152263.3', gene = gen1, exon_start = 1, exon_start_genomic = ex1,
                                  exon_end = 8, exon_end_genomic = ex2)
        tr2 = Transcript_region(transcript = 'NM_002529.3', gene = gen2, exon_start = 10, exon_start_genomic= ex3,
                                  exon_end = 17, exon_end_genomic = ex4)
        molecular_fusion = Fusion(r_frame_preserved = True, domains = domain, components = [tr1, linker_seq, tr2], causative_event = cause_eve)
        print(molecular_fusion.json())
    except ValidationError as f:
        print(f)

    filename = 'example_gene_fusion_class.json'
    with open(filename, 'w') as files:
        files.write(molecular_fusion.make_json())
    files.close()


if __name__ == '__main__':
    main()

