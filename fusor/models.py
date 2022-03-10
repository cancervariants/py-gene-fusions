"""Model for fusion class"""
from typing import Optional, List, Union, Literal
from enum import Enum
from abc import ABC

from pydantic import BaseModel, validator, StrictInt, StrictBool, StrictStr, \
    Extra, ValidationError, root_validator
from ga4gh.vrsatile.pydantic import return_value
from ga4gh.vrsatile.pydantic.vrsatile_models import GeneDescriptor, \
    LocationDescriptor, SequenceDescriptor, CURIE
from ga4gh.vrsatile.pydantic.vrs_models import Sequence, SequenceLocation


class BaseModelForbidExtra(BaseModel):
    """Base model with extra fields forbidden."""

    class Config:
        """Configure class."""

        extra = Extra.forbid


class FUSORTypes(str, Enum):
    """Define FUSOR object type values."""

    FUNCTIONAL_DOMAIN = "FunctionalDomain"
    TRANSCRIPT_SEGMENT_COMPONENT = "TranscriptSegmentComponent"
    TEMPLATED_SEQUENCE_COMPONENT = "TemplatedSequenceComponent"
    LINKER_SEQUENCE_COMPONENT = "LinkerSequenceComponent"
    GENE_COMPONENT = "GeneComponent"
    UNKNOWN_GENE_COMPONENT = "UnknownGeneComponent"
    ANY_GENE_COMPONENT = "AnyGeneComponent"
    REGULATORY_ELEMENT = "RegulatoryElement"
    FUSION = "Fusion"


class AdditionalFields(str, Enum):
    """Define possible fields that can be added to Fusion object."""

    SEQUENCE_ID = "sequence_id"
    LOCATION_ID = "location_id"
    GENE_DESCRIPTOR = "gene_descriptor"


class DomainStatus(str, Enum):
    """Define possible statuses of functional domains."""

    LOST = "lost"
    PRESERVED = "preserved"


class FunctionalDomain(BaseModel):
    """Define FunctionalDomain class"""

    type: Literal[FUSORTypes.FUNCTIONAL_DOMAIN] = FUSORTypes.FUNCTIONAL_DOMAIN
    id: CURIE
    name: StrictStr
    status: DomainStatus
    gene_descriptor: GeneDescriptor
    location_descriptor: LocationDescriptor

    _get_id_val = validator("id", allow_reuse=True)(return_value)

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "FunctionalDomain",
                "status": "lost",
                "name": "Tyrosine-protein kinase, catalytic domain",
                "id": "interpro:IPR020635",
                "gene_descriptor": {
                    "id": "gene:NTRK1",
                    "gene_id": "hgnc:8031",
                    "label": "8031",
                    "type": "GeneDescriptor",
                },
                "location_descriptor": {
                    "id": "fusor.location_descriptor:NP_002520.2",
                    "type": "LocationDescriptor",
                    "location": {
                        "sequence_id": "ga4gh:SQ.vJvm06Wl5J7DXHynR9ksW7IK3_3jlFK6",  # noqa: E501
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
            }


class Component(ABC, BaseModel):
    """Define base component class."""

    pass


class TranscriptSegmentComponent(Component):
    """Define TranscriptSegment class"""

    type: Literal[FUSORTypes.TRANSCRIPT_SEGMENT_COMPONENT] = FUSORTypes.TRANSCRIPT_SEGMENT_COMPONENT  # noqa: E501
    transcript: CURIE
    exon_start: Optional[StrictInt]
    exon_start_offset: Optional[StrictInt] = 0
    exon_end: Optional[StrictInt]
    exon_end_offset: Optional[StrictInt] = 0
    gene_descriptor: GeneDescriptor
    component_genomic_start: Optional[LocationDescriptor]
    component_genomic_end: Optional[LocationDescriptor]

    @root_validator(pre=True)
    def check_exons(cls, values):
        """Check that at least one of {`exon_start`, `exon_end`} is set.
        If set, check that the corresponding `component_genomic` field is set.
        If not set, set corresponding offset to `None`

        """
        msg = "Must give values for either `exon_start`, `exon_end`, or both"
        exon_start = values.get("exon_start")
        exon_end = values.get("exon_end")
        assert exon_start or exon_end, msg

        if exon_start:
            msg = "Must give `component_genomic_start` if `exon_start` is given"  # noqa: E501
            assert values.get("component_genomic_start"), msg
        else:
            values["exon_start_offset"] = None

        if exon_end:
            msg = "Must give `component_genomic_end` if `exon_end` is given"
            assert values.get("component_genomic_end"), msg
        else:
            values["exon_end_offset"] = None
        return values

    _get_transcript_val = validator("transcript", allow_reuse=True)(return_value)  # noqa: E501

    def nomenclature(self, i: Optional[int] = None) -> str:
        """Generate component nomenclature."""
        prefix = f"{self.transcript}({self.gene_descriptor.label})"
        start, start_offset, end, end_offset = "", "", "", ""
        if i != -1:
            end = self.exon_end
            if self.exon_end_offset is not None:
                end_offset = self.exon_end_offset
        if i != 0:
            start = self.exon_start
            if self.exon_start_offset is not None:
                start_offset = self.exon_start_offset
        return f"{prefix}: e.{start}{start_offset}_{end}{end_offset}"

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "TranscriptSegmentComponent",
                "transcript": "refseq:NM_152263.3",
                "exon_start": 1,
                "exon_start_offset": 0,
                "exon_end": 8,
                "exon_end_offset": 0,
                "gene_descriptor": {
                    "type": "GeneDescriptor",
                    "id": "gene:TPM3",
                    "gene_id": "hgnc:12012",
                    "type": "GeneDescriptor",
                    "label": "TPM3",
                },
                "component_genomic_start": {
                    "id": "TPM3:exon1",
                    "type": "LocationDescriptor",
                    "location_id": "ga4gh:VSL.vyyyExx4enSZdWZr3z67-T8uVKH50uLi",  # noqa: E501
                    "location": {
                        "sequence_id": "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT",  # noqa: E501
                        "type": "SequenceLocation",
                        "interval": {
                            "start": {
                                "type": "Number",
                                "value": 154192135
                            },
                            "end": {
                                "type": "Number",
                                "value": 154192136
                            },
                            "type": "SequenceInterval"
                        }
                    }
                },
                "component_genomic_end": {
                    "id": "TPM3:exon8",
                    "type": "LocationDescriptor",
                    "location_id": "ga4gh:VSL._1bRdL4I6EtpBvVK5RUaXb0NN3k0gpqa",  # noqa: E501
                    "location": {
                        "sequence_id": "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT",  # noqa: E501
                        "type": "SequenceLocation",
                        "interval": {
                            "start": {
                                "type": "Number",
                                "value": 154170398
                            },
                            "end": {
                                "type": "Number",
                                "value": 154170399
                            },
                            "type": "SequenceInterval"
                        }
                    }
                }
            }


class LinkerComponent(Component):
    """Define Linker class (linker sequence)"""

    type: Literal[FUSORTypes.LINKER_SEQUENCE_COMPONENT] = FUSORTypes.LINKER_SEQUENCE_COMPONENT  # noqa: E501
    linker_sequence: SequenceDescriptor

    def nomenclature(self, _: Optional[int] = None) -> str:
        """Generate component nomenclature."""
        return self.linker_sequence.sequence  # type: ignore

    @validator("linker_sequence", pre=True)
    def validate_sequence(cls, v):
        """Enforce nucleotide base code requirements on sequence literals."""
        if isinstance(v, dict):
            try:
                v["sequence"] = v["sequence"].upper()
                seq = v["sequence"]
            except KeyError:
                raise TypeError
        elif isinstance(v, SequenceDescriptor):
            v.sequence = v.sequence.upper()
            seq = v.sequence
        else:
            raise TypeError

        try:
            Sequence(__root__=seq)
        except ValidationError:
            raise AssertionError("sequence does not match regex '^[A-Za-z*\\-]*$'")  # noqa: E501

        return v

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "LinkerSequenceComponent",
                "linker_sequence": {
                    "id": "sequence:ACGT",
                    "type": "SequenceDescriptor",
                    "sequence": "ACGT",
                    "residue_type": "SO:0000348"
                }
            }


class Strand(str, Enum):
    """Define possible values for strand"""

    POSITIVE = "+"
    NEGATIVE = "-"


class TemplatedSequenceComponent(BaseModel):
    """Define Templated Sequence Component class.
    A templated sequence is a contiguous genomic sequence found in the
    gene product
    """

    type: Literal[FUSORTypes.TEMPLATED_SEQUENCE_COMPONENT] = FUSORTypes.TEMPLATED_SEQUENCE_COMPONENT  # noqa: E501
    region: LocationDescriptor
    strand: Strand

    # add strand to sequencelocation, add chr property
    def nomenclature(self, _: Optional[int] = None) -> str:
        """Generate component nomenclature."""
        if self.region and self.region.location:
            location = self.region.location
            if isinstance(location, SequenceLocation):
                sequence_id = location.sequence_id
                start = location.interval.start.number
                end = location.interval.end.number
                return f"{sequence_id}: g.{start}_{end}({self.strand.value})"
            else:
                raise ValueError
        else:
            raise ValueError

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "TemplatedSequenceComponent",
                "region": {
                    "id": "chr12:44908821-44908822(+)",
                    "type": "LocationDescriptor",
                    "location_id": "ga4gh:VSL.AG54ZRBhg6pwpPLafF4KgaAHpdFio6l5",  # noqa: E501
                    "location": {
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl",  # noqa: E501
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 44908821},
                            "end": {"type": "Number", "value": 44908822}
                        },
                    },
                    "label": "chr12:44908821-44908822(+)"
                },
                "strand": "+"
            }


class GeneComponent(Component):
    """Define Gene component class."""

    type: Literal[FUSORTypes.GENE_COMPONENT] = FUSORTypes.GENE_COMPONENT
    gene_descriptor: GeneDescriptor

    def nomenclature(self, _: Optional[int] = None) -> str:
        """Generate component nomenclature."""
        if self.gene_descriptor.gene_id:
            gene_id = self.gene_descriptor.gene_id
        elif self.gene_descriptor.gene and self.gene_descriptor.gene.gene_id:
            gene_id = self.gene_descriptor.gene.gene_id
        else:
            raise ValueError
        return f"{self.gene_descriptor.label}({gene_id})"

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "GeneComponent",
                "gene_descriptor": {
                    "id": "gene:BRAF",
                    "gene_id": "hgnc:1097",
                    "label": "BRAF",
                    "type": "GeneDescriptor",
                }
            }


class UnknownGeneComponent(BaseModel):
    """Define UnknownGene class. This is primarily intended to represent a
    partner in the result of a fusion partner-agnostic assay, which identifies
    the absence of an expected gene. For example, a FISH break-apart probe may
    indicate rearrangement of an MLL gene, but by design, the test cannot
    provide the identity of the new partner. In this case, we would associate
    any clinical observations from this patient with the fusion of MLL with
    an UnknownGene component.
    """

    type: Literal[FUSORTypes.UNKNOWN_GENE_COMPONENT] = FUSORTypes.UNKNOWN_GENE_COMPONENT  # noqa: E501

    def nomenclature(self, _: Optional[int] = None) -> str:
        """Generate component nomenclature."""
        return "?"

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "UnknownGeneComponent"
            }


class AnyGeneComponent(BaseModel):
    """Define AnyGene class. This is primarily intended to represent a partner
    in a categorical fusion, typifying generalizable characteristics of a class
    of fusions such as retained or lost regulatory elements and/or functional
    domains, often curated from biomedical literature for use in genomic
    knowledgebases. For example, EWSR1 rearrangements are often found in Ewing
    and Ewing-like small round cell sarcomas, regardless of the partner gene.
    We would associate this assertion with the fusion of EWSR1 with an
    AnyGene component.
    """

    type: Literal[FUSORTypes.ANY_GENE_COMPONENT] = FUSORTypes.ANY_GENE_COMPONENT  # noqa: E501

    def nomenclature(self, _: Optional[int] = None) -> str:
        """Generate component nomenclature."""
        return "*"

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "AnyGeneComponent"
            }


class Event(str, Enum):
    """Define Event class (causative event)"""

    REARRANGEMENT = "rearrangement"
    READTHROUGH = "read-through"
    TRANSSPLICING = "trans-splicing"


class RegulatoryElementType(str, Enum):
    """Define possible types of Regulatory Elements. Options are the possible
    values for /regulatory_class value property in the INSDC controlled
    vocabulary.
    https://www.insdc.org/controlled-vocabulary-regulatoryclass
    """

    ATTENUATOR = "attenuator"
    CAAT_SIGNAL = "caat_signal"
    ENHANCER = "enhancer"
    ENHANCER_BLOCKING_ELEMENT = "enhancer_blocking_element"
    GC_SIGNAL = "gc_signal"
    IMPRINTING_CONTROL_REGION = "imprinting_control_region"
    INSULATOR = "insulator"
    LOCUS_CONTROL_REGION = "locus_control_region"
    MINUS_35_SIGNAL = "minus_35_signal"
    MINUS_10_SIGNAL = "minus_10_signal"
    POLYA_SIGNAL_SEQUENCE = "polya_signal_sequence"
    PROMOTER = "promoter"
    RESPONSE_ELEMENT = "response_element"
    RIBOSOME_BINDING_SITE = "ribosome_binding_site"
    RIBOSWITCH = "riboswitch"
    SILENCER = "silencer"
    TATA_BOX = "tata_box"
    TERMINATOR = "terminator"
    OTHER = "other"


class RegulatoryElement(BaseModel):
    """Define RegulatoryElement class"""

    type: Literal[FUSORTypes.REGULATORY_ELEMENT] = FUSORTypes.REGULATORY_ELEMENT  # noqa: E501
    element_type: RegulatoryElementType
    reference_id: Optional[CURIE] = None
    gene_descriptor: Optional[GeneDescriptor] = None
    genomic_location_start: Optional[LocationDescriptor] = None
    genomic_lcoation_end: Optional[LocationDescriptor] = None

    _get_ref_id_val = validator("reference_id", allow_reuse=True)(return_value)

    @root_validator(pre=True)
    def ensure_min_values(cls, values):
        """Ensure one of {`reference_id`, `gene_descriptor`,
        `location_descriptor`} is set.
        """
        if not any([values["reference_id"], values["gene_descriptor"],
                    values["genomic_location_start"] and values["genomic_location_end"]]):  # noqa: E501
            raise ValueError
        return values

    def nomenclature(self, i: Optional[int]) -> str:
        """Generate component nomenclature representation.
        :param Optional[int] i: index position in fusion component list. If
            1st, should == 0. If last, should == -1. Leave unset otherwise.
        :return: string with nomenclature representation of component.
        """
        nm_string = f"reg_{self.type.value}"
        if self.reference_id:
            nm_string += f"_{self.reference_id}"
        # if element reference:
        # _<reference_id>
        # if located element:
        # _<Chromosome ID>(chr <1-22,X,Y>): 
        # g.<start coordinate>_<end coordinate>

        nm_string += f"@{self.gene_descriptor.gene_id}"

        return nm_string

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "RegulatoryElement",
                "element_type": "Promoter",
                "gene_descriptor": {
                    "id": "gene:BRAF",
                    "gene_id": "hgnc:1097",
                    "label": "BRAF",
                    "type": "GeneDescriptor",
                }
            }


class Fusion(BaseModel, ABC):
    """Define Fusion class"""

    regulatory_elements: Optional[List[RegulatoryElement]]
    structural_components: List[Component]

    @validator("structural_components")
    def structural_components_length(cls, v):
        """Ensure >=2 structural components"""
        if len(v) < 2:
            raise ValueError("Fusion must contain at least 2 structural "
                             "components.")
        else:
            return v

    @validator("structural_components")
    def structural_components_ends(cls, v):
        """Ensure start/end components are of legal types and have fields
        required by their position.
        """
        # TODO first check for regulatoryelement
        if isinstance(v[0], TranscriptSegmentComponent):
            if v[0].exon_start is None:
                raise ValueError
        elif isinstance(v[0], LinkerComponent):
            raise ValueError

        if isinstance(v[-1], TranscriptSegmentComponent):
            if v[-1].exon_end is None:
                raise ValueError
        elif isinstance(v[-1], LinkerComponent):
            raise ValueError


class Evidence(str, Enum):
    """Form of evidence supporting identification of the fusion."""

    OBSERVED = "observed"
    INFERRED = "inferred"


class MolecularAssay(BaseModelForbidExtra):
    """Information pertaining to the assay used in identifying the fusion."""

    assay_name: Optional[StrictStr]
    eco_id: Optional[CURIE]
    eco_label: Optional[StrictStr]
    method_uri: Optional[CURIE]

    _get_eco_id_val = validator("eco_id", allow_reuse=True)(return_value)
    _get_method_id_val = validator("method_id", allow_reuse=True)(return_value)

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            schema["example"] = {
                "method_uri": "pmid:34974290",
                "eco_id": "ECO:0005629",
                "eco_label": "DNA affinity chromatography evidence used in manual assertion",  # noqa: E501
                "assay_name": "Breakapart FISH probe"
            }


class AssayedFusion(Fusion):
    """Assayed gene fusions from biological specimens are directly detected
    using RNA-based gene fusion assays, or alternatively may be inferred from
    genomic rearrangements detected by whole genome sequencing or by
    coarser-scale cytogenomic assays. Example: an EWSR1 fusion inferred
    from a breakapart FISH assay.
    """

    causative_event: Optional[Event]
    fusion_evidence: Optional[Evidence]
    molecular_assay: Optional[MolecularAssay]
    structural_components: List[Union[TranscriptSegmentComponent,
                                      GeneComponent,
                                      TemplatedSequenceComponent,
                                      LinkerComponent,
                                      UnknownGeneComponent]]

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "causative_event": "rearrangement",
                "fusion_evidence": "observed",
                "molecular_assay": '{"ECO_id": "ECO:0007159"}',
                "structural_components": [
                    {
                        "component_type": "gene",
                        "gene_descriptor": {
                            "id": "gene:EWSR1",
                            "gene_id": "hgnc:3058",
                            "label": "EWSR1",
                            "type": "GeneDescriptor",
                        }
                    },
                    {
                        "component_type": "UnknownGene"
                    }
                ]
            }


class CategoricalFusion(Fusion):
    """Categorical gene fusions are generalized concepts representing a class
    of fusions by their shared attributes, such as retained or lost regulatory
    elements and/or functional domains, and are typically curated from the
    biomedical literature for use in genomic knowledgebases.
    """

    r_frame_preserved: Optional[StrictBool]
    functional_domains: Optional[List[FunctionalDomain]]
    structural_components: List[Union[TranscriptSegmentComponent,
                                      GeneComponent,
                                      TemplatedSequenceComponent,
                                      LinkerComponent,
                                      AnyGeneComponent]]

    class Config(BaseModelForbidExtra.Config):
        """Configure class."""

        @staticmethod
        def schema_extra(schema, _):
            """Provide example"""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "type": "Fusion",
                "r_frame_preserved": True,
                "functional_domains": [
                    {
                        "status": "lost",
                        "name": "cystatin domain",
                        "id": "interpro:IPR000010",
                        "gene": {
                            "id": "gene:CST1",
                            "gene_id": "hgnc:2743",
                            "label": "CST1",
                            "type": "GeneDescriptor",
                        }
                    }
                ],
                "structural_components": [
                    {
                        "component_type": "TranscriptSegment",
                        "transcript": "refseq:NM_152263.3",
                        "exon_start": 1,
                        "exon_start_offset": 0,
                        "exon_end": 8,
                        "exon_end_offset": 0,
                        "gene": {
                            "id": "gene:TPM3",
                            "gene_id": "hgnc:12012",
                            "type": "GeneDescriptor",
                            "label": "TPM3",
                        },
                        "component_genomic_start": {
                            "id": "TPM3:exon1",
                            "type": "LocationDescriptor",
                            "location_id": "ga4gh:VSL.vyyyExx4enSZdWZr3z67-T8uVKH50uLi",  # noqa: E501
                            "location": {
                                "sequence_id": "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT",  # noqa: E501
                                "type": "SequenceLocation",
                                "interval": {
                                    "start": {
                                        "type": "Number",
                                        "value": 154192135
                                    },
                                    "end": {
                                        "type": "Number",
                                        "value": 154192136
                                    },
                                    "type": "SequenceInterval"
                                }
                            }
                        },
                        "component_genomic_end": {
                            "id": "TPM3:exon8",
                            "type": "LocationDescriptor",
                            "location_id": "ga4gh:VSL._1bRdL4I6EtpBvVK5RUaXb0NN3k0gpqa",  # noqa: E501
                            "location": {
                                "sequence_id": "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT",  # noqa: E501
                                "type": "SequenceLocation",
                                "interval": {
                                    "start": {
                                        "type": "Number",
                                        "value": 154170398
                                    },
                                    "end": {
                                        "type": "Number",
                                        "value": 154170399
                                    },
                                    "type": "SequenceInterval"
                                }
                            }
                        }
                    },
                    {
                        "component_type": "gene",
                        "gene": {
                            "id": "gene:ALK",
                            "type": "GeneDescriptor",
                            "gene_id": "hgnc:427",
                            "label": "ALK"
                        }
                    }
                ],
                "regulatory_elements": [
                    {
                        "element_type": "promoter",
                        "gene": {
                            "id": "gene:BRAF",
                            "type": "GeneDescriptor",
                            "gene_id": "hgnc:1097",
                            "label": "BRAF"
                        }
                    }
                ]
            }
