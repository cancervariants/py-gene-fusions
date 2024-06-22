"""Model for fusion class"""
from abc import ABC
from enum import Enum
from typing import Any, Literal

from ga4gh.vrsatile.pydantic import return_value
from ga4gh.vrsatile.pydantic.vrsatile_models import (
    CURIE,
    GeneDescriptor,
    LocationDescriptor,
    SequenceDescriptor,
)
from pydantic import (
    BaseModel,
    ConfigDict,
    StrictBool,
    StrictInt,
    StrictStr,
    field_validator,
    model_validator,
)
from pydantic.fields import Field


class BaseModelForbidExtra(BaseModel, extra="forbid"):
    """Base Pydantic model class with extra values forbidden."""


class FUSORTypes(str, Enum):
    """Define FUSOR object type values."""

    FUNCTIONAL_DOMAIN = "FunctionalDomain"
    TRANSCRIPT_SEGMENT_ELEMENT = "TranscriptSegmentElement"
    TEMPLATED_SEQUENCE_ELEMENT = "TemplatedSequenceElement"
    LINKER_SEQUENCE_ELEMENT = "LinkerSequenceElement"
    GENE_ELEMENT = "GeneElement"
    UNKNOWN_GENE_ELEMENT = "UnknownGeneElement"
    MULTIPLE_POSSIBLE_GENES_ELEMENT = "MultiplePossibleGenesElement"
    REGULATORY_ELEMENT = "RegulatoryElement"
    CATEGORICAL_FUSION = "CategoricalFusion"
    ASSAYED_FUSION = "AssayedFusion"
    CAUSATIVE_EVENT = "CausativeEvent"


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
    status: DomainStatus
    associated_gene: GeneDescriptor
    id: CURIE | None = Field(None, alias="_id")
    label: StrictStr | None = None
    sequence_location: LocationDescriptor | None = None

    _get_id_val = field_validator("id")(return_value)

    model_config = ConfigDict(
        populate_by_name=True,
        json_schema_extra={
            "example": {
                "type": "FunctionalDomain",
                "status": "lost",
                "label": "Tyrosine-protein kinase, catalytic domain",
                "_id": "interpro:IPR020635",
                "associated_gene": {
                    "id": "gene:NTRK1",
                    "gene_id": "hgnc:8031",
                    "label": "8031",
                    "type": "GeneDescriptor",
                },
                "sequence_location": {
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
            }
        },
    )


class StructuralElementType(str, Enum):
    """Define possible structural element type values."""

    TRANSCRIPT_SEGMENT_ELEMENT = FUSORTypes.TRANSCRIPT_SEGMENT_ELEMENT.value
    TEMPLATED_SEQUENCE_ELEMENT = FUSORTypes.TEMPLATED_SEQUENCE_ELEMENT.value
    LINKER_SEQUENCE_ELEMENT = FUSORTypes.LINKER_SEQUENCE_ELEMENT.value
    GENE_ELEMENT = FUSORTypes.GENE_ELEMENT.value
    UNKNOWN_GENE_ELEMENT = FUSORTypes.UNKNOWN_GENE_ELEMENT.value
    MULTIPLE_POSSIBLE_GENES_ELEMENT = FUSORTypes.MULTIPLE_POSSIBLE_GENES_ELEMENT.value


class BaseStructuralElement(ABC, BaseModel):
    """Define base structural element class."""

    type: StructuralElementType


class TranscriptSegmentElement(BaseStructuralElement):
    """Define TranscriptSegment class"""

    type: Literal[
        FUSORTypes.TRANSCRIPT_SEGMENT_ELEMENT
    ] = FUSORTypes.TRANSCRIPT_SEGMENT_ELEMENT
    transcript: CURIE
    exon_start: StrictInt | None = None
    exon_start_offset: StrictInt | None = 0
    exon_end: StrictInt | None = None
    exon_end_offset: StrictInt | None = 0
    gene_descriptor: GeneDescriptor
    element_genomic_start: LocationDescriptor | None = None
    element_genomic_end: LocationDescriptor | None = None

    @model_validator(mode="before")
    def check_exons(cls, values):
        """Check that at least one of {`exon_start`, `exon_end`} is set.
        If set, check that the corresponding `element_genomic` field is set.
        If not set, set corresponding offset to `None`

        """
        msg = "Must give values for either `exon_start`, `exon_end`, or both"
        exon_start = values.get("exon_start")
        exon_end = values.get("exon_end")
        if (not exon_start) and (not exon_end):
            raise ValueError(msg)

        if exon_start:
            if not values.get("element_genomic_start"):
                msg = "Must give `element_genomic_start` if `exon_start` is given"
                raise ValueError(msg)
        else:
            values["exon_start_offset"] = None

        if exon_end:
            if not values.get("element_genomic_end"):
                msg = "Must give `element_genomic_end` if `exon_end` is given"
                raise ValueError(msg)
        else:
            values["exon_end_offset"] = None
        return values

    _get_transcript_val = field_validator("transcript")(return_value)
    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_152263.3",
                "exon_start": 1,
                "exon_start_offset": 0,
                "exon_end": 8,
                "exon_end_offset": 0,
                "gene_descriptor": {
                    "id": "normalize.gene:TPM3",
                    "type": "GeneDescriptor",
                    "label": "TPM3",
                    "gene_id": "hgnc:12012",
                },
                "element_genomic_start": {
                    "id": "fusor.location_descriptor:NC_000001.11",
                    "type": "LocationDescriptor",
                    "label": "NC_000001.11",
                    "location": {
                        "type": "SequenceLocation",
                        "sequence_id": "refseq:NC_000001.11",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 154192135},
                            "end": {"type": "Number", "value": 154192136},
                        },
                    },
                },
                "element_genomic_end": {
                    "id": "fusor.location_descriptor:NC_000001.11",
                    "type": "LocationDescriptor",
                    "label": "NC_000001.11",
                    "location": {
                        "type": "SequenceLocation",
                        "sequence_id": "refseq:NC_000001.11",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 154170399},
                            "end": {"type": "Number", "value": 154170400},
                        },
                    },
                },
            }
        },
    )


class LinkerElement(BaseStructuralElement, extra="forbid"):
    """Define Linker class (linker sequence)"""

    type: Literal[
        FUSORTypes.LINKER_SEQUENCE_ELEMENT
    ] = FUSORTypes.LINKER_SEQUENCE_ELEMENT
    linker_sequence: SequenceDescriptor

    @field_validator("linker_sequence", mode="before")
    def validate_sequence(cls, v):
        """Enforce nucleotide base code requirements on sequence literals."""
        if isinstance(v, dict):
            try:
                v["sequence"] = v["sequence"].upper()
            except KeyError as e:
                raise TypeError from e
        elif isinstance(v, SequenceDescriptor):
            v.sequence = v.sequence.upper()
        else:
            raise TypeError

        return v

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "type": "LinkerSequenceElement",
                "linker_sequence": {
                    "id": "sequence:ACGT",
                    "type": "SequenceDescriptor",
                    "sequence": "ACGT",
                    "residue_type": "SO:0000348",
                },
            }
        },
    )


class Strand(str, Enum):
    """Define possible values for strand"""

    POSITIVE = "+"
    NEGATIVE = "-"


class TemplatedSequenceElement(BaseStructuralElement):
    """Define Templated Sequence Element class.
    A templated sequence is a contiguous genomic sequence found in the gene
    product.
    """

    type: Literal[
        FUSORTypes.TEMPLATED_SEQUENCE_ELEMENT
    ] = FUSORTypes.TEMPLATED_SEQUENCE_ELEMENT
    region: LocationDescriptor
    strand: Strand

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "type": "TemplatedSequenceElement",
                "region": {
                    "id": "chr12:44908821-44908822(+)",
                    "type": "LocationDescriptor",
                    "location_id": "ga4gh:VSL.AG54ZRBhg6pwpPLafF4KgaAHpdFio6l5",
                    "location": {
                        "type": "SequenceLocation",
                        "sequence_id": "ga4gh:SQ.6wlJpONE3oNb4D69ULmEXhqyDZ4vwNfl",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 44908821},
                            "end": {"type": "Number", "value": 44908822},
                        },
                    },
                    "label": "chr12:44908821-44908822(+)",
                },
                "strand": "+",
            }
        },
    )


class GeneElement(BaseStructuralElement):
    """Define Gene Element class."""

    type: Literal[FUSORTypes.GENE_ELEMENT] = FUSORTypes.GENE_ELEMENT
    gene_descriptor: GeneDescriptor

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "type": "GeneElement",
                "gene_descriptor": {
                    "id": "gene:BRAF",
                    "gene_id": "hgnc:1097",
                    "label": "BRAF",
                    "type": "GeneDescriptor",
                },
            }
        },
    )


class UnknownGeneElement(BaseStructuralElement):
    """Define UnknownGene class. This is primarily intended to represent a
    partner in the result of a fusion partner-agnostic assay, which identifies
    the absence of an expected gene. For example, a FISH break-apart probe may
    indicate rearrangement of an MLL gene, but by design, the test cannot
    provide the identity of the new partner. In this case, we would associate
    any clinical observations from this patient with the fusion of MLL with
    an UnknownGene element.
    """

    type: Literal[FUSORTypes.UNKNOWN_GENE_ELEMENT] = FUSORTypes.UNKNOWN_GENE_ELEMENT

    model_config = ConfigDict(
        json_schema_extra={"example": {"type": "UnknownGeneElement"}},
    )


class MultiplePossibleGenesElement(BaseStructuralElement):
    """Define MultiplePossibleGenesElement class. This is primarily intended to
    represent a partner in a categorical fusion, typifying generalizable
    characteristics of a class of fusions such as retained or lost regulatory elements
    and/or functional domains, often curated from biomedical literature for use in
    genomic knowledgebases. For example, EWSR1 rearrangements are often found in Ewing
    and Ewing-like small round cell sarcomas, regardless of the partner gene.
    We would associate this assertion with the fusion of EWSR1 with a
    MultiplePossibleGenesElement.
    """

    type: Literal[
        FUSORTypes.MULTIPLE_POSSIBLE_GENES_ELEMENT
    ] = FUSORTypes.MULTIPLE_POSSIBLE_GENES_ELEMENT

    model_config = ConfigDict(
        json_schema_extra={"example": {"type": "MultiplePossibleGenesElement"}},
    )


class RegulatoryClass(str, Enum):
    """Define possible classes of Regulatory Elements. Options are the possible values
    for /regulatory_class value property in the INSDC controlled vocabulary:
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
    """Define RegulatoryElement class.

    `feature_id` would ideally be constrained as a CURIE, but Encode, our preferred
    feature ID source, doesn't currently have a registered CURIE structure for EH_
    identifiers. Consequently, we permit any kind of free text.
    """

    type: Literal[FUSORTypes.REGULATORY_ELEMENT] = FUSORTypes.REGULATORY_ELEMENT
    regulatory_class: RegulatoryClass
    feature_id: str | None = None
    associated_gene: GeneDescriptor | None = None
    feature_location: LocationDescriptor | None = None

    _get_ref_id_val = field_validator("feature_id")(return_value)

    @model_validator(mode="before")
    def ensure_min_values(cls, values):
        """Ensure that one of {`feature_id`, `feature_location`}, and/or
        `associated_gene` is set.
        """
        if not (
            bool(values.get("feature_id")) ^ bool(values.get("feature_location"))
        ) and not (values.get("associated_gene")):
            msg = "Must set 1 of {`feature_id`, `associated_gene`} and/or `feature_location`"
            raise ValueError(msg)
        return values

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "type": "RegulatoryElement",
                "regulatory_class": "promoter",
                "feature_location": {
                    "type": "LocationDescriptor",
                    "id": "fusor.location_descriptor:NC_000001.11",
                    "location": {
                        "sequence_id": "ga4gh:SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        "type": "SequenceLocation",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 155593},
                            "end": {"type": "Number", "value": 155610},
                        },
                    },
                },
            }
        },
    )


class FusionType(str, Enum):
    """Specify possible Fusion types."""

    CATEGORICAL_FUSION = FUSORTypes.CATEGORICAL_FUSION.value
    ASSAYED_FUSION = FUSORTypes.ASSAYED_FUSION.value

    @classmethod
    def values(cls) -> set:
        """Provide all possible enum values."""
        return {c.value for c in cls}


class AbstractFusion(BaseModel, ABC):
    """Define Fusion class"""

    type: FusionType
    regulatory_element: RegulatoryElement | None = None
    structural_elements: list[BaseStructuralElement]

    @classmethod
    def _access_object_attr(
        cls,
        obj: dict | BaseModel,
        attr_name: str,
    ) -> Any | None:  # noqa: ANN401
        """Help enable safe access of object properties while performing
        validation for Pydantic class objects. Because the validator could be handling either
        existing Pydantic class objects, or candidate dictionaries, we need a flexible
        accessor function.

        :param obj: object to access
        :param attr_name: name of attribute to retrieve
        :return: attribute if successful retrieval, otherwise None
        :raise ValueError: if object doesn't have properties (ie it's not a dict or Pydantic
            model)
        """
        if isinstance(obj, BaseModel):
            try:
                return obj.__getattribute__(attr_name)
            except AttributeError:
                return None
        elif isinstance(obj, dict):
            return obj.get(attr_name)
        else:
            msg = "Unrecognized type, should only pass entities with properties"
            raise ValueError(msg)

    @classmethod
    def _fetch_gene_id(
        cls,
        obj: dict | BaseModel,
        gene_descriptor_field: str,
    ) -> str | None:
        """Get gene ID if element includes a gene annotation.

        :param obj: element to fetch gene from. Might not contain a gene (e.g. it's a
            TemplatedSequenceElement) so we have to use safe checks to fetch.
        :param gene_descriptor_field: name of gene_descriptor field
        :return: gene ID if gene is defined
        """
        gene_descriptor = cls._access_object_attr(obj, gene_descriptor_field)
        if gene_descriptor:
            gene_value = cls._access_object_attr(gene_descriptor, "gene")
            if gene_value:
                gene_id = cls._access_object_attr(gene_value, "gene_id")
                if gene_id:
                    return gene_id
            gene_id = cls._access_object_attr(gene_descriptor, "gene_id")
            if gene_id:
                return gene_id
        return None

    @model_validator(mode="before")
    def enforce_abc(cls, values):
        """Ensure only subclasses can be instantiated."""
        if cls.__name__ == "AbstractFusion":
            msg = "Cannot instantiate Fusion abstract class"
            raise ValueError(msg)
        return values

    @model_validator(mode="before")
    def enforce_element_quantities(cls, values):
        """Ensure minimum # of elements, and require > 1 unique genes.

        To validate the unique genes rule, we extract gene IDs from the elements that
        designate genes, and take the number of total elements. If there is only one
        unique gene ID, and there are no non-gene-defining elements (such as
        an unknown partner), then we raise an error.
        """
        qt_error_msg = (
            "Fusions must contain >= 2 structural elements, or >=1 structural element "
            "and a regulatory element"
        )
        structural_elements = values.get("structural_elements", [])
        if not structural_elements:
            raise ValueError(qt_error_msg)
        num_structural_elements = len(structural_elements)
        reg_element = values.get("regulatory_element")
        if (num_structural_elements + bool(reg_element)) < 2:
            raise ValueError(qt_error_msg)

        uq_gene_msg = "Fusions must form a chimeric transcript from two or more genes, or a novel interaction between a rearranged regulatory element with the expressed product of a partner gene."
        gene_ids = []
        if reg_element:
            gene_id = cls._fetch_gene_id(
                obj=reg_element, gene_descriptor_field="associated_gene"
            )
            if gene_id:
                gene_ids.append(gene_id)

        for element in structural_elements:
            gene_id = cls._fetch_gene_id(
                obj=element, gene_descriptor_field="gene_descriptor"
            )
            if gene_id:
                gene_ids.append(gene_id)

        unique_gene_ids = set(gene_ids)
        if len(unique_gene_ids) == 1 and len(gene_ids) == (
            num_structural_elements + bool(reg_element)
        ):
            raise ValueError(uq_gene_msg)
        return values

    @model_validator(mode="after")
    def structural_elements_ends(cls, values):
        """Ensure start/end elements are of legal types and have fields
        required by their position.
        """
        elements = values.structural_elements
        if isinstance(elements[0], TranscriptSegmentElement):
            if elements[0].exon_end is None and not values["regulatory_element"]:
                msg = "5' TranscriptSegmentElement fusion partner must contain ending exon position"
                raise ValueError(msg)
        elif isinstance(elements[0], LinkerElement):
            msg = "First structural element cannot be LinkerSequence"
            raise ValueError(msg)

        if len(elements) > 2:
            for element in elements[1:-1]:
                if isinstance(element, TranscriptSegmentElement) and (
                    element.exon_start is None or element.exon_end is None
                ):
                    msg = "Connective TranscriptSegmentElement must include both start and end positions"
                    raise ValueError(msg)
        if isinstance(elements[-1], TranscriptSegmentElement) and (
            elements[-1].exon_start is None
        ):
            msg = "3' fusion partner junction must include " "starting position"
            raise ValueError
        return values


class Evidence(str, Enum):
    """Form of evidence supporting identification of the fusion."""

    OBSERVED = "observed"
    INFERRED = "inferred"


class Assay(BaseModelForbidExtra):
    """Information pertaining to the assay used in identifying the fusion."""

    type: Literal["Assay"] = "Assay"
    assay_name: StrictStr | None = None
    assay_id: CURIE | None = None
    method_uri: CURIE | None = None
    fusion_detection: Evidence | None = None

    _get_assay_id_val = field_validator("assay_id")(return_value)
    _get_method_uri_val = field_validator("method_uri")(return_value)

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "method_uri": "pmid:33576979",
                "assay_id": "obi:OBI_0003094",
                "assay_name": "fluorescence in-situ hybridization assay",
                "fusion_detection": "inferred",
            }
        }
    )


AssayedFusionElements = list[
    TranscriptSegmentElement
    | GeneElement
    | TemplatedSequenceElement
    | LinkerElement
    | UnknownGeneElement
]


class EventType(str, Enum):
    """Permissible values for describing the underlying causative event driving an
    assayed fusion.
    """

    REARRANGEMENT = "rearrangement"
    READ_THROUGH = "read-through"
    TRANS_SPLICING = "trans-splicing"


class CausativeEvent(BaseModelForbidExtra):
    """The evaluation of a fusion may be influenced by the underlying mechanism that
    generated the fusion. Often this will be a DNA rearrangement, but it could also be
    a read-through or trans-splicing event.
    """

    type: Literal[FUSORTypes.CAUSATIVE_EVENT] = FUSORTypes.CAUSATIVE_EVENT
    event_type: EventType
    event_description: StrictStr | None = None

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "type": "CausativeEvent",
                "event_type": "rearrangement",
                "event_description": "chr2:g.pter_8,247,756::chr11:g.15,825,273_cen_qter (der11) and chr11:g.pter_15,825,272::chr2:g.8,247,757_cen_qter (der2)",
            }
        },
    )


class AssayedFusion(AbstractFusion):
    """Assayed gene fusions from biological specimens are directly detected using
    RNA-based gene fusion assays, or alternatively may be inferred from genomic
    rearrangements detected by whole genome sequencing or by coarser-scale cytogenomic
    assays. Example: an EWSR1 fusion inferred from a breakapart FISH assay.
    """

    type: Literal[FUSORTypes.ASSAYED_FUSION] = FUSORTypes.ASSAYED_FUSION
    structural_elements: AssayedFusionElements
    causative_event: CausativeEvent | None = None
    assay: Assay | None = None

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "type": "AssayedFusion",
                "causative_event": {
                    "type": "CausativeEvent",
                    "event_type": "rearrangement",
                    "event_description": "chr2:g.pter_8,247,756::chr11:g.15,825,273_cen_qter (der11) and chr11:g.pter_15,825,272::chr2:g.8,247,757_cen_qter (der2)",
                },
                "assay": {
                    "type": "Assay",
                    "method_uri": "pmid:33576979",
                    "assay_id": "obi:OBI_0003094",
                    "assay_name": "fluorescence in-situ hybridization assay",
                    "fusion_detection": "inferred",
                },
                "structural_elements": [
                    {
                        "type": "GeneElement",
                        "gene_descriptor": {
                            "id": "gene:EWSR1",
                            "gene_id": "hgnc:3058",
                            "label": "EWSR1",
                            "type": "GeneDescriptor",
                        },
                    },
                    {"type": "UnknownGeneElement"},
                ],
            }
        },
    )


CategoricalFusionElements = list[
    TranscriptSegmentElement
    | GeneElement
    | TemplatedSequenceElement
    | LinkerElement
    | MultiplePossibleGenesElement
]


class CategoricalFusion(AbstractFusion):
    """Categorical gene fusions are generalized concepts representing a class
    of fusions by their shared attributes, such as retained or lost regulatory
    elements and/or functional domains, and are typically curated from the
    biomedical literature for use in genomic knowledgebases.
    """

    type: Literal[FUSORTypes.CATEGORICAL_FUSION] = FUSORTypes.CATEGORICAL_FUSION
    r_frame_preserved: StrictBool | None = None
    critical_functional_domains: list[FunctionalDomain] | None = None
    structural_elements: CategoricalFusionElements

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "type": "CategoricalFusion",
                "r_frame_preserved": True,
                "critical_functional_domains": [
                    {
                        "type": "FunctionalDomain",
                        "status": "lost",
                        "label": "cystatin domain",
                        "id": "interpro:IPR000010",
                        "associated_gene": {
                            "id": "gene:CST1",
                            "gene_id": "hgnc:2743",
                            "label": "CST1",
                            "type": "GeneDescriptor",
                        },
                    }
                ],
                "structural_elements": [
                    {
                        "type": "TranscriptSegmentElement",
                        "transcript": "refseq:NM_152263.3",
                        "exon_start": 1,
                        "exon_start_offset": 0,
                        "exon_end": 8,
                        "exon_end_offset": 0,
                        "gene_descriptor": {
                            "id": "gene:TPM3",
                            "gene_id": "hgnc:12012",
                            "type": "GeneDescriptor",
                            "label": "TPM3",
                        },
                        "element_genomic_start": {
                            "id": "TPM3:exon1",
                            "type": "LocationDescriptor",
                            "location_id": "ga4gh:VSL.vyyyExx4enSZdWZr3z67-T8uVKH50uLi",
                            "location": {
                                "sequence_id": "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT",
                                "type": "SequenceLocation",
                                "interval": {
                                    "start": {"type": "Number", "value": 154192135},
                                    "end": {"type": "Number", "value": 154192136},
                                    "type": "SequenceInterval",
                                },
                            },
                        },
                        "element_genomic_end": {
                            "id": "TPM3:exon8",
                            "type": "LocationDescriptor",
                            "location_id": "ga4gh:VSL._1bRdL4I6EtpBvVK5RUaXb0NN3k0gpqa",
                            "location": {
                                "sequence_id": "ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT",
                                "type": "SequenceLocation",
                                "interval": {
                                    "start": {"type": "Number", "value": 154170398},
                                    "end": {"type": "Number", "value": 154170399},
                                    "type": "SequenceInterval",
                                },
                            },
                        },
                    },
                    {
                        "type": "GeneElement",
                        "gene_descriptor": {
                            "id": "gene:ALK",
                            "type": "GeneDescriptor",
                            "gene_id": "hgnc:427",
                            "label": "ALK",
                        },
                    },
                ],
                "regulatory_element": {
                    "type": "RegulatoryElement",
                    "regulatory_class": "promoter",
                    "associated_gene": {
                        "id": "gene:BRAF",
                        "type": "GeneDescriptor",
                        "gene_id": "hgnc:1097",
                        "label": "BRAF",
                    },
                },
            }
        },
    )


Fusion = CategoricalFusion | AssayedFusion
