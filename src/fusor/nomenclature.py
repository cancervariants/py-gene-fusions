"""Provide helper methods for fusion nomenclature generation."""

from biocommons.seqrepo.seqrepo import SeqRepo
from ga4gh.vrs.models import SequenceLocation

from fusor.exceptions import IDTranslationException
from fusor.models import (
    GeneElement,
    RegulatoryClass,
    RegulatoryElement,
    TemplatedSequenceElement,
    TranscriptSegmentElement,
)
from fusor.tools import translate_identifier


def reg_element_nomenclature(element: RegulatoryElement, sr: SeqRepo) -> str:
    """Return fusion nomenclature for regulatory element.

    :param element: a regulatory element object
    :param sr: a SeqRepo instance
    :return: regulatory element nomenclature representation
    :raises ValueError: if unable to retrieve genomic location or coordinates,
        or if missing element reference ID, genomic location, and associated
        gene
    """
    element_class = element.regulatoryClass.value
    if element_class == RegulatoryClass.ENHANCER:
        type_string = "e"
    elif element_class == RegulatoryClass.PROMOTER:
        type_string = "p"
    else:
        type_string = f"{element.regulatoryClass.value}"
    feature_string = ""
    if element.featureId:
        feature_string += f"_{element.featureId}"
    elif element.featureLocation:
        feature_location = element.featureLocation
        sequence_id = feature_location.sequenceReference.id
        refseq_id = str(translate_identifier(sr, sequence_id, "refseq")).split(":")[1]
        try:
            chrom = str(translate_identifier(sr, sequence_id, "GRCh38")).split(":")[1]
        except IDTranslationException as e:
            raise ValueError from e
        feature_string += f"_{refseq_id}(chr {chrom}):g.{feature_location.start}_{feature_location.end}"
    if element.associatedGene:
        if element.associatedGene.id:
            gene_id = element.associatedGene.id
        else:
            raise ValueError
        feature_string += f"@{element.associatedGene.label}({gene_id})"
    if not feature_string:
        raise ValueError
    return f"reg_{type_string}{feature_string}"


def tx_segment_nomenclature(element: TranscriptSegmentElement) -> str:
    """Return fusion nomenclature for transcript segment element

    :param element: a tx segment element. Treated as a junction component if only one
        end is provided.
    :return: element nomenclature representation
    """
    transcript = str(element.transcript)
    if ":" in transcript:
        transcript = transcript.split(":")[1]

    prefix = f"{transcript}({element.gene.label})"
    start = element.exonStart if element.exonStart else ""
    if element.exonStartOffset:
        if element.exonStartOffset > 0:
            start_offset = f"+{element.exonStartOffset}"
        else:
            start_offset = str(element.exonStartOffset)
    else:
        start_offset = ""
    end = element.exonEnd if element.exonEnd else ""
    if element.exonEndOffset:
        if element.exonEndOffset > 0:
            end_offset = f"+{element.exonEndOffset}"
        else:
            end_offset = str(element.exonEndOffset)
    else:
        end_offset = ""
    return f"{prefix}:e.{start}{start_offset}{'_' if start and end else ''}{end}{end_offset}"


def templated_seq_nomenclature(element: TemplatedSequenceElement, sr: SeqRepo) -> str:
    """Return fusion nomenclature for templated sequence element.

    :param element: a templated sequence element
    :param sr: SeqRepo instance to use
    :return: element nomenclature representation
    :raises ValueError: if location isn't a SequenceLocation or if unable
        to retrieve region or location
    """
    region = element.region
    if region:
        sequence_reference = element.region.sequenceReference
        if isinstance(sequence_reference, SequenceLocation):
            sequence_id = str(sequence_reference.id)
            refseq_id = str(translate_identifier(sr, sequence_id, "refseq"))
            start = region.start
            end = region.end
            try:
                chrom = str(translate_identifier(sr, sequence_id, "GRCh38")).split(":")[
                    1
                ]
            except IDTranslationException as e:
                raise ValueError from e
            return f"{refseq_id.split(':')[1]}(chr {chrom}):g.{start}_{end}({element.strand.value})"
        raise ValueError
    raise ValueError


def gene_nomenclature(element: GeneElement) -> str:
    """Return fusion nomenclature for gene element.

    :param element: a gene element object
    :return: element nomenclature representation
    :raises ValueError: if unable to retrieve gene ID
    """
    if element.gene.id:
        gene_id = element.gene.id
    else:
        raise ValueError
    return f"{element.gene.label}({gene_id})"
