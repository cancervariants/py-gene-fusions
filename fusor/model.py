"""Model for fusion class"""
import json
from pydantic import BaseModel, validator, StrictInt, StrictBool, StrictStr, \
    Extra
from typing import Optional, List, Union, Literal
from enum import Enum
from ga4gh.vrsatile.pydantic import return_value
from ga4gh.vrsatile.pydantic.vrsatile_model import GeneDescriptor, \
    LocationDescriptor, SequenceDescriptor, CURIE


class DomainStatus(str, Enum):
    """Define possible statuses of critical domains."""

    LOST = "lost"
    PRESERVED = "preserved"


class CriticalDomain(BaseModel):
    """Define CriticalDomain class"""

    status: DomainStatus
    name: StrictStr
    id: CURIE
    gene_descriptor: GeneDescriptor

    _get_id_val = validator('id', allow_reuse=True)(return_value)

    class Config:
        """Configure class."""

        extra = Extra.forbid

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                'status': 'lost',
                'name': 'cystatin domain',
                'id': 'interpro:IPR000010',
                'gene_descriptor': {
                    'id': 'gene:CST1',
                    'gene_id': 'hgnc:2743',
                    'label': 'CST1',
                    'type': 'GeneDescriptor',
                }
            }


class ComponentType(str, Enum):
    """Define possible transcript components."""

    TRANSCRIPT_SEGMENT = 'transcript_segment'
    GENOMIC_REGION = 'genomic_region'
    LINKER_SEQUENCE = 'linker_sequence'
    GENE = 'gene'
    UNKNOWN_GENE = 'unknown_gene'


class TranscriptSegmentComponent(BaseModel):
    """Define TranscriptSegment class"""

    component_type: Literal[ComponentType.TRANSCRIPT_SEGMENT] = ComponentType.TRANSCRIPT_SEGMENT  # noqa: E501
    transcript: CURIE
    exon_start: StrictInt
    exon_start_offset: StrictInt = 0
    exon_end: StrictInt
    exon_end_offset: StrictInt = 0
    gene_descriptor: GeneDescriptor
    component_genomic_region: LocationDescriptor

    _get_transcript_val = validator('transcript', allow_reuse=True)(return_value)  # noqa: E501

    class Config:
        """Configure class."""

        extra = Extra.forbid

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                'component_type': 'transcript_segment',
                'transcript': 'refseq:NM_152263.3',
                'exon_start': 1,
                'exon_start_offset': 0,
                'exon_end': 8,
                'exon_end_offset': 0,
                'gene_descriptor': {
                    'id': 'gene:TPM3',
                    'gene_id': 'hgnc:12012',
                    'type': 'GeneDescriptor',
                    'label': 'TPM3',
                },
                'component_genomic_region': {
                    'id': 'TPM3:exon1-exon8',
                    'type': 'LocationDescriptor',
                    'location': {
                        'sequence_id': 'ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT',  # noqa: E501
                        'type': 'SequenceLocation',
                        'interval': {
                            'start': {
                                'type': 'Number',
                                'value': 154192135,
                            },
                            'end': {
                                'type': 'Number',
                                'value': 154170399,
                            },
                            'type': 'SequenceInterval',
                        }
                    }
                }
            }


class LinkerComponent(BaseModel):
    """Define Linker class (linker sequence)"""

    component_type: Literal[ComponentType.LINKER_SEQUENCE] = ComponentType.LINKER_SEQUENCE  # noqa: E501
    linker_sequence: SequenceDescriptor

    @validator('linker_sequence')
    def validate_sequence(cls, v):
        """Enforce nucleotide base code requirements on sequence literals."""
        if isinstance(v, dict):
            try:
                sequence = v['sequence']
            except KeyError:
                raise TypeError
        elif isinstance(v, SequenceDescriptor):
            sequence = v.sequence
        else:
            raise TypeError
        msg = 'Linker sequence must consist only of {A,C,G,T}'
        assert set('ACGT') >= set(sequence), msg
        return v

    class Config:
        """Configure class."""

        extra = Extra.forbid

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                'component_type': 'linker_sequence',
                'linker_sequence': {
                    'id': 'sequence:ACGT',
                    'type': 'SequenceDescriptor',
                    'sequence': 'ACGT',
                    'residue_type': 'SO:0000348'
                }
            }


class Strand(str, Enum):
    """Define possible values for strand"""

    POSITIVE = "+"
    NEGATIVE = "-"


class GenomicRegionComponent(BaseModel):
    """Define GenomicRegion component class."""

    component_type: Literal[ComponentType.GENOMIC_REGION] = ComponentType.GENOMIC_REGION  # noqa: E501
    region: LocationDescriptor
    strand: Strand

    # add strand to sequencelocation, add chr property

    class Config:
        """Configure class."""

        extra = Extra.forbid

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                'component_type': 'genomic_region',
                'region': {
                    'id': 'chr12:44908821-44908822(+)',
                    'type': 'LocationDescriptor',
                    'location': {
                        'type': 'SequenceLocation',
                        'sequence_id': 'ga4gh:SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl',  # noqa: E501
                        'interval': {
                            'type': 'SequenceInterval',
                            'start': {'type': 'Number', 'value': 44908821},
                            'end': {'type': 'Number', 'value': 44908822},
                        },
                    },
                    'label': 'chr12:44908821-44908822(+)'
                },
                'strand': '+'
            }


class GeneComponent(BaseModel):
    """Define Gene component class."""

    component_type: Literal[ComponentType.GENE] = ComponentType.GENE
    gene_descriptor: GeneDescriptor

    class Config:
        """Configure class."""

        extra = Extra.forbid

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                'component_type': 'gene',
                'gene_descriptor': {
                    'id': 'gene:BRAF',
                    'gene_id': 'hgnc:1097',
                    'label': 'BRAF',
                    'type': 'GeneDescriptor',
                }
            }


class UnknownGeneComponent(BaseModel):
    """Define UnknownGene class"""

    component_type: Literal[ComponentType.UNKNOWN_GENE] = ComponentType.UNKNOWN_GENE  # noqa: E501

    class Config:
        """Configure class."""

        extra = Extra.forbid

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                'component_type': 'unknown_gene'
            }


class Event(str, Enum):
    """Define Event class (causative event)"""

    REARRANGEMENT = 'rearrangement'
    READTHROUGH = 'read-through'
    TRANSSPLICING = 'trans-splicing'


class RegulatoryElementType(str, Enum):
    """Define possible types of Regulatory Elements."""

    PROMOTER = 'promoter'
    ENHANCER = 'enhancer'


class RegulatoryElement(BaseModel):
    """Define RegulatoryElement class"""

    type: RegulatoryElementType
    gene_descriptor: GeneDescriptor

    class Config:
        """Configure class."""

        extra = Extra.forbid

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                'type': 'promoter',
                'gene_descriptor': {
                    'id': 'gene:BRAF',
                    'gene_id': 'hgnc:1097',
                    'label': 'BRAF',
                    'type': 'GeneDescriptor',
                }
            }


class Fusion(BaseModel):
    """Define Fusion class"""

    r_frame_preserved: Optional[StrictBool]
    protein_domains: Optional[List[CriticalDomain]]
    transcript_components: List[Union[TranscriptSegmentComponent,
                                      GeneComponent,
                                      GenomicRegionComponent,
                                      LinkerComponent,
                                      UnknownGeneComponent]]
    causative_event: Optional[Event]
    regulatory_elements: Optional[List[RegulatoryElement]]

    @validator('transcript_components')
    def transcript_components_length(cls, v):
        """Ensure >=2 transcript components"""
        if len(v) < 2:
            raise ValueError('Fusion must contain at least 2 transcript '
                             'components.')
        else:
            return v

    def make_json(self):
        """JSON helper function"""
        return json.dumps(self, default=lambda o: o.__dict__, indent=4)

    class Config:
        """Configure class."""

        extra = Extra.forbid

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                'r_frame_preserved': True,
                'protein_domains': [
                    {
                        'status': 'lost',
                        'name': 'cystatin domain',
                        'id': 'interpro:IPR000010',
                        'gene': {
                            'id': 'gene:CST1',
                            'gene_id': 'hgnc:2743',
                            'label': 'CST1',
                            'type': 'GeneDescriptor',
                        }
                    }
                ],
                'transcript_components': [
                    {
                        'component_type': 'transcript_segment',
                        'transcript': 'refseq:NM_152263.3',
                        'exon_start': 1,
                        'exon_start_offset': 0,
                        'exon_end': 8,
                        'exon_end_offset': 0,
                        'gene': {
                            'id': 'gene:TPM3',
                            'gene_id': 'hgnc:12012',
                            'type': 'GeneDescriptor',
                            'label': 'TPM3',
                        },
                        'component_genomic_region': {
                            'id': 'TPM3:exon1-exon8',
                            'type': 'LocationDescriptor',
                            'location': {
                                'sequence_id': 'ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT',  # noqa: E501
                                'type': 'SequenceLocation',
                                'interval': {
                                    'start': {
                                        'type': 'Number',
                                        'value': 154192135
                                    },
                                    'end': {
                                        'type': 'Number',
                                        'value': 154170399
                                    },
                                    'type': 'SequenceInterval'
                                }
                            }
                        }
                    },
                    {
                        'component_type': 'gene',
                        'gene': {
                            'id': 'gene:ALK',
                            'type': 'GeneDescriptor',
                            'gene_id': 'hgnc:427',
                            'label': 'ALK'
                        }
                    }
                ],
                'causative_event': 'rearrangement',
                'regulatory_elements': [
                    {
                        'type': 'promoter',
                        'gene': {
                            'id': 'gene:BRAF',
                            'type': 'GeneDescriptor',
                            'gene_id': 'hgnc:1097',
                            'label': 'BRAF'
                        }
                    }
                ]
            }
