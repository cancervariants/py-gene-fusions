"""Model for fusion class"""

from abc import ABC
from enum import Enum
from typing import Any, Literal

from ga4gh.core.domain_models import Gene
from ga4gh.vrs.models import (
    LiteralSequenceExpression,
    SequenceLocation,
)
from gene.schemas import CURIE
from pydantic import (
    BaseModel,
    ConfigDict,
    StrictBool,
    StrictInt,
    StrictStr,
    model_validator,
)


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


class DomainStatus(str, Enum):
    """Define possible statuses of functional domains."""

    LOST = "lost"
    PRESERVED = "preserved"


class FunctionalDomain(BaseModel):
    """Define FunctionalDomain class"""

    type: Literal[FUSORTypes.FUNCTIONAL_DOMAIN] = FUSORTypes.FUNCTIONAL_DOMAIN
    status: DomainStatus
    gene: Gene
    id: CURIE | None
    label: StrictStr | None = None
    sequenceLocation: SequenceLocation | None = None

    model_config = ConfigDict(
        populate_by_name=True,
        json_schema_extra={
            "example": {
                "type": "FunctionalDomain",
                "status": "lost",
                "label": "Tyrosine-protein kinase, catalytic domain",
                "id": "interpro:IPR020635",
                "gene": {
                    "id": "gene:NTRK1",
                    "gene_id": "hgnc:8031",
                    "label": "8031",
                    "type": "Gene",
                },
                "sequenceLocation": {
                    "id": "SL.vJvm06Wl5J7DXHynR9ksW7IK3_3jlFK6",
                    "start": 510,
                    "end": "781",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "id": "NP_002520.2",
                        "label": "GRCh38:chr22",
                        "type": "SequenceReference",
                        "refgetAccession": "SQ.7B7SHsmchAR0dFcDCuSFjJAo7tX87krQ",
                        "residueAlphabet": "na",
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

    type: Literal[FUSORTypes.TRANSCRIPT_SEGMENT_ELEMENT] = (
        FUSORTypes.TRANSCRIPT_SEGMENT_ELEMENT
    )
    transcript: CURIE
    exonStart: StrictInt | None = None
    exonStartOffset: StrictInt | None = 0
    exonEnd: StrictInt | None = None
    exonEndOffset: StrictInt | None = 0
    gene: Gene
    elementGenomicStart: SequenceLocation | None = None
    elementGenomicEnd: SequenceLocation | None = None

    @model_validator(mode="before")
    def check_exons(cls, values):
        """Check that at least one of {`exonStart`, `exonEnd`} is set.
        If set, check that the corresponding `elementGenomic` field is set.
        If not set, set corresponding offset to `None`

        """
        msg = "Must give values for either `exonStart`, `exonEnd`, or both"
        exon_start = values.get("exonStart")
        exon_end = values.get("exonEnd")
        if (not exon_start) and (not exon_end):
            raise ValueError(msg)

        if exon_start:
            if not values.get("elementGenomicStart"):
                msg = "Must give `elementGenomicStart` if `exonStart` is given"
                raise ValueError(msg)
        else:
            values["exonStartOffset"] = None

        if exon_end:
            if not values.get("elementGenomicEnd"):
                msg = "Must give `elementGenomicEnd` if `exonEnd` is given"
                raise ValueError(msg)
        else:
            values["exonEndOffset"] = None
        return values

    model_config = ConfigDict(
        # TODO: verify this example
        json_schema_extra={
            "example": {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_152263.3",
                "exonStart": 1,
                "exonStartOffset": 0,
                "exonEnd": 8,
                "exonEndOffset": 0,
                "gene": {
                    "id": "hgnc:12012",
                    "type": "Gene",
                    "label": "TPM3",
                },
                "elementGenomicStart": {
                    # TODO: this digest may not be correct, but putting a placeholder example here for now
                    "id": "SL.vJvm06Wl5J7DXHynR9ksW7IK3_3jlFK6",
                    "type": "SequenceLocation",
                    "label": "NC_000001.11",
                    "sequenceReference": {
                        "type": "SequenceReference",
                        "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                    },
                    "start": 154192135,
                    "end": 154192136,
                },
                "elementGenomicEnd": {
                    # TODO: this digest may not be correct, but putting a placeholder example here for now
                    "id": "SL.vJvm06Wl5J7DXHynR9ksW7IK3_3jlFK6",
                    "type": "SequenceLocation",
                    "label": "NC_000001.11",
                    "sequenceReference": {
                        "type": "SequenceReference",
                        "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                    },
                    "start": 154170399,
                    "end": 154170400,
                },
            }
        },
    )


class LinkerElement(BaseStructuralElement, extra="forbid"):
    """Define Linker class (linker sequence)"""

    type: Literal[FUSORTypes.LINKER_SEQUENCE_ELEMENT] = (
        FUSORTypes.LINKER_SEQUENCE_ELEMENT
    )
    linkerSequence: LiteralSequenceExpression

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        # TODO: update this example json once models approved
        json_schema_extra={
            "example": {
                "type": "LinkerSequenceElement",
                "linkerSequence": {
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

    type: Literal[FUSORTypes.TEMPLATED_SEQUENCE_ELEMENT] = (
        FUSORTypes.TEMPLATED_SEQUENCE_ELEMENT
    )
    region: SequenceLocation
    strand: Strand

    model_config = ConfigDict(
        # TODO: update this example json once models approved
        json_schema_extra={
            "example": {
                "type": "TemplatedSequenceElement",
                "region": {
                    "id": "chr12:44908821-44908822(+)",
                    "type": "SequenceLocation",
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
    gene: Gene

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "type": "GeneElement",
                "gene": {
                    "id": "hgnc:1097",
                    "label": "BRAF",
                    "type": "Gene",
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

    type: Literal[FUSORTypes.MULTIPLE_POSSIBLE_GENES_ELEMENT] = (
        FUSORTypes.MULTIPLE_POSSIBLE_GENES_ELEMENT
    )

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

    `featureId` would ideally be constrained as a CURIE, but Encode, our preferred
    feature ID source, doesn't currently have a registered CURIE structure for EH_
    identifiers. Consequently, we permit any kind of free text.
    """

    type: Literal[FUSORTypes.REGULATORY_ELEMENT] = FUSORTypes.REGULATORY_ELEMENT
    regulatoryClass: RegulatoryClass
    featureId: str | None = None
    associatedGene: Gene | None = None
    featureLocation: SequenceLocation | None = None

    @model_validator(mode="before")
    def ensure_min_values(cls, values):
        """Ensure that one of {`featureId`, `featureLocation`}, and/or
        `associatedGene` is set.
        """
        if not (
            bool(values.get("featureId")) ^ bool(values.get("featureLocation"))
        ) and not (values.get("associatedGene")):
            msg = (
                "Must set 1 of {`featureId`, `associatedGene`} and/or `featureLocation`"
            )
            raise ValueError(msg)
        return values

    model_config = ConfigDict(
        # TODO: update this example json once models approved
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
    regulatoryElement: RegulatoryElement | None = None
    structure: list[BaseStructuralElement]
    readingFramePreserved: StrictBool | None = None

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
    ) -> str | None:
        """Get gene ID if element includes a gene annotation.

        :param obj: element to fetch gene from. Might not contain a gene (e.g. it's a
            TemplatedSequenceElement) so we have to use safe checks to fetch.
        :return: gene ID if gene is defined
        """
        gene_info = cls._access_object_attr(obj, "gene")
        if gene_info:
            gene_id = cls._access_object_attr(gene_info, "id")
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
        structure = values.get("structure", [])
        if not structure:
            raise ValueError(qt_error_msg)
        num_structure = len(structure)
        reg_element = values.get("regulatoryElement")
        if (num_structure + bool(reg_element)) < 2:
            raise ValueError(qt_error_msg)

        uq_gene_msg = "Fusions must form a chimeric transcript from two or more genes, or a novel interaction between a rearranged regulatory element with the expressed product of a partner gene."
        gene_ids = []
        if reg_element:
            gene_id = cls._fetch_gene_id(obj=reg_element)
            if gene_id:
                gene_ids.append(gene_id)

        for element in structure:
            gene_id = cls._fetch_gene_id(obj=element)
            if gene_id:
                gene_ids.append(gene_id)

        unique_gene_ids = set(gene_ids)
        if len(unique_gene_ids) == 1 and len(gene_ids) == (
            num_structure + bool(reg_element)
        ):
            raise ValueError(uq_gene_msg)
        return values

    @model_validator(mode="after")
    def structure_ends(cls, values):
        """Ensure start/end elements are of legal types and have fields
        required by their position.
        """
        elements = values.structure
        if isinstance(elements[0], TranscriptSegmentElement):
            if elements[0].exonEnd is None and not values["regulatory_element"]:
                msg = "5' TranscriptSegmentElement fusion partner must contain ending exon position"
                raise ValueError(msg)
        elif isinstance(elements[0], LinkerElement):
            msg = "First structural element cannot be LinkerSequence"
            raise ValueError(msg)

        if len(elements) > 2:
            for element in elements[1:-1]:
                if isinstance(element, TranscriptSegmentElement) and (
                    element.exonStart is None or element.exonEnd is None
                ):
                    msg = "Connective TranscriptSegmentElement must include both start and end positions"
                    raise ValueError(msg)
        if isinstance(elements[-1], TranscriptSegmentElement) and (
            elements[-1].exonStart is None
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
    assayName: StrictStr | None = None
    assayId: CURIE | None = None
    methodUri: CURIE | None = None
    fusionDetection: Evidence | None = None

    model_config = ConfigDict(
        # TODO: verify this example json once models approved
        json_schema_extra={
            "example": {
                "methodUri": "pmid:33576979",
                "assayId": "obi:OBI_0003094",
                "assayName": "fluorescence in-situ hybridization assay",
                "fusionDetection": "inferred",
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
    eventType: EventType
    eventDescription: StrictStr | None = None

    model_config = ConfigDict(
        # TODO: verify this example json once models approved
        json_schema_extra={
            "example": {
                "type": "CausativeEvent",
                "eventType": "rearrangement",
                "eventDescription": "chr2:g.pter_8,247,756::chr11:g.15,825,273_cen_qter (der11) and chr11:g.pter_15,825,272::chr2:g.8,247,757_cen_qter (der2)",
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
    structure: AssayedFusionElements
    causativeEvent: CausativeEvent | None = None
    assay: Assay | None = None

    model_config = ConfigDict(
        # TODO: update this example json once models approved
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
                "structure": [
                    {
                        "type": "GeneElement",
                        "gene": {
                            "type": "Gene",
                            "id": "gene:EWSR1",
                            "gene_id": "hgnc:3058",
                            "label": "EWSR1",
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
    criticalFunctionalDomains: list[FunctionalDomain] | None = None
    structure: CategoricalFusionElements

    model_config = ConfigDict(
        # TODO: update this example json once models approved
        json_schema_extra={
            "example": {
                "type": "CategoricalFusion",
                "reading_frame_preserved": True,
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
                "structure": [
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
                        "type": "Gene",
                        "id": "gene:ALK",
                        "gene_id": "hgnc:427",
                        "label": "ALK",
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
