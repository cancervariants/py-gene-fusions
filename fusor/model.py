"""Model for fusion class"""
import json
from pydantic import BaseModel, validator, StrictInt, StrictBool
from typing import Optional, List, Union
from gene.schemas import GeneDescriptor, SequenceLocation, ChromosomeLocation
from enum import Enum


def check_curie(cls, v):
    """Validate curies."""
    if v is not None:
        def _is_curie(value):
            """Check that value is a curie

            :param str value: Value to validate
            """
            assert all(
                [
                    value.count(':') == 1,
                    value.find(' ') == -1,
                    value[-1] != ':'
                ]
            ), 'must be a CURIE'

        if isinstance(v, str):
            _is_curie(v)
        elif isinstance(v, list):
            for item in v:
                _is_curie(item)
    return v


class DomainStatus(str, Enum):
    """Define possible statuses of critical domains."""

    LOST = "lost"
    PRESERVED = "preserved"


class CriticalDomain(BaseModel):
    """Define CriticalDomain class"""

    status: DomainStatus
    name: str
    id: str
    gene: GeneDescriptor

    _validate_id = validator('id', allow_reuse=True)(check_curie)

    class Config:
        """Configure class."""

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
                'gene': {
                    'id': 'gene:CST1',
                    'value_id': 'hgnc:2743',
                    'label': 'CST1'
                }
            }


class TranscriptComponent(BaseModel):
    """Define abstract TranscriptComponentClass."""

    component_type: str


class GenomicRegion(BaseModel):
    """Define GenomicRegion class"""

    type = 'LocationDescription'
    description: Optional[str] = None
    value: Union[SequenceLocation, ChromosomeLocation]


class TranscriptSegmentComponent(TranscriptComponent):
    """Define TranscriptSegment class"""

    component_type = 'transcript_segment'
    transcript: str
    exon_start: StrictInt
    exon_start_offset: StrictInt = 0
    exon_end: StrictInt
    exon_end_offset: StrictInt = 0
    gene: GeneDescriptor
    component_genomic_region: GenomicRegion

    _validate_transcript = \
        validator('transcript', allow_reuse=True)(check_curie)

    class Config:
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                'component_type': 'transcript_segment',
                'transcript': 'NM_152263.3',
                'exon_start': 1,
                'exon_start_offset': 0,
                'exon_end': 8,
                'exon_end_offset': 0,
                'gene': 'TPM3',
                'component_genomic_region': {
                    'type': 'LocationDescription',
                    'value': {
                        'sequence_id': '',
                        'type': 'SequenceLocation',
                        'interval': {
                            'start': 154192135,
                            'end': 154170399,
                            'type': 'SimpleInterval'
                        }
                    }
                }
            }


class LinkerComponent(TranscriptComponent):
    """Define Linker class (linker sequence)"""

    component_type = 'linker_sequence'
    sequence: str

    @validator('sequence')
    def valid_input(cls, v):
        """Validate linker_sequence"""
        assert set(v.upper()) <= set('ACGT'), 'Linker sequence must ' \
                                              'only contain A,C,G,T'
        return v

    class Config:
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                'component_type': 'linker_sequence',
                'sequence': 'AGCTC'
            }


class Strand(Enum):
    """Define valid Strand values."""

    POSITIVE = '+'
    NEGATIVE = '-'


class GenomicRegionComponent(TranscriptComponent):
    """Define GenomicRegion component class."""

    component_type = 'genomic_region'
    chr: str
    strand: Strand
    start: StrictInt
    end: StrictInt

    class Config:
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                'component_type': 'genomic_region',
                'chr': '12',
                'strand': '+',
                'start': 3498,
                'end': 3500
            }


class GeneComponent(TranscriptComponent):
    """Define Gene component class."""

    component_type = 'gene'
    gene: GeneDescriptor

    class Config:
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                'component_type': 'gene',
                'gene': {
                    'id': 'gene:BRAF',
                    'value_id': 'hgnc:1097',
                    'label': 'BRAF'
                }
            }


class UnknownGeneComponent(TranscriptComponent):
    """Define UnknownGene class"""

    component_type = 'unknown_gene'

    class Config:
        """Configure class."""

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


class Event(Enum):
    """Define Event class (causative event)"""

    REARRANGEMENT = 'rearrangement'
    READTHROUGH = 'read-through'
    TRANSSPLICING = 'trans-splicing'


class RegulatoryElementType(Enum):
    """Define possible types of Regulatory Elements."""

    PROMOTER = 'promoter'
    ENHANCER = 'enhancer'


class RegulatoryElement(BaseModel):
    """Define RegulatoryElement class"""

    type: RegulatoryElementType
    gene: GeneDescriptor

    class Config:
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if 'title' in schema.keys():
                schema.pop('title', None)
            for prop in schema.get('properties', {}).values():
                prop.pop('title', None)
            schema['example'] = {
                'type': 'promoter',
                'gene': {
                    'id': 'gene:BRAF',
                    'value_id': 'hgnc:1097',
                    'label': 'BRAF'
                }
            }


class Fusion(BaseModel):
    """Define Fusion class"""

    r_frame_preserved: Optional[StrictBool]
    protein_domains: Optional[List[CriticalDomain]]
    transcript_components: List[TranscriptComponent]
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
                            'value_id': 'hgnc:2743',
                            'label': 'CST1'
                        }
                    }
                ],
                'transcript_components': [
                    {
                        'component_type': 'transcript_segment',
                        'transcript': 'NM_152263.3',
                        'exon_start': 1,
                        'exon_start_offset': 0,
                        'exon_end': 8,
                        'exon_end_offset': 0,
                        'gene': 'TPM3',
                        'component_genomic_region': {
                            'type': 'LocationDescription',
                            'value': {
                                'sequence_id': '',
                                'type': 'SequenceLocation',
                                'interval': {
                                    'start': 154192135,
                                    'end': 154170399,
                                    'type': 'SimpleInterval'
                                }
                            }
                        }
                    },
                    {
                        'component_type': 'gene',
                        'gene': {
                            'id': 'gene:ALK',
                            'value_id': 'hgnc:427',
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
                            'value_id': 'hgnc:1097',
                            'label': 'BRAF'
                        }
                    }
                ]
            }
