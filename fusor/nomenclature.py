"""Provide fusion nomenclature generation methods."""
from biocommons.seqrepo.seqrepo import SeqRepo
from ga4gh.vrsatile.pydantic.vrs_models import SequenceLocation
from fusor.exceptions import IDTranslationException

from fusor.models import GeneComponent, RegulatoryElement, \
    TemplatedSequenceComponent, TranscriptSegmentComponent
from fusor.tools import translate_identifier


def reg_element_nomenclature(element: RegulatoryElement, sr: SeqRepo) -> str:
    """Return fusion nomenclature for regulatory element.
    :param RegulatoryElement element: a regulatory element object
    :param SeqRepo sr: a SeqRepo instance
    :return: regulatory element nomenclature representation
    :raises ValueError: if unable to retrieve genomic location or coordinates,
        or if missing element reference ID, genomic location, and associated
        gene
    """
    nm_type_string = f"reg_{element.type.value}"
    nm_string = ""
    if element.element_reference:
        nm_string += f"_{element.element_reference}"
    elif element.genomic_location:
        start = element.genomic_location
        sequence_id = start.location.sequence_id
        try:
            chr = str(
                translate_identifier(sr, sequence_id, "GRCh38")
            ).split(":")[1]
        except IDTranslationException:
            raise ValueError
        nm_string += f"_{sequence_id}(chr{chr}):g.{start.location.interval.start.value}_{start.location.interval.end.value}"  # noqa: E501
    if element.associated_gene:
        if element.associated_gene.gene_id:
            gene_id = gene_id = element.associated_gene.gene_id

        if element.associated_gene.gene_id:
            gene_id = element.associated_gene.gene_id
        elif element.associated_gene.gene and \
                element.associated_gene.gene.gene_id:
            gene_id = element.associated_gene.gene.gene_id
        else:
            raise ValueError
        nm_string += f"@{element.associated_gene.label}({gene_id})"
    if not nm_string:
        raise ValueError
    return nm_type_string + nm_string


def tx_segment_nomenclature(component: TranscriptSegmentComponent,
                            i: int) -> str:
    """Return fusion nomenclature for transcript segment component
    :param TranscriptSegmentComponent component: a tx segment component
    :param int i: index of component in fusion components sequence
    :return: component nomenclature representation
    """
    prefix = f"{component.transcript}({component.gene_descriptor.label})"
    start, start_offset, end, end_offset = "", "", "", ""
    if i != -1:
        end = component.exon_end
        if component.exon_end_offset is not None:
            end_offset = component.exon_end_offset
    if i != 0:
        start = component.exon_start
        if component.exon_start_offset is not None:
            start_offset = component.exon_start_offset
    return f"{prefix}: e.{start}{start_offset}_{end}{end_offset}"


def templated_seq_nomenclature(component: TemplatedSequenceComponent) -> str:
    """Return fusion nomenclature for templated sequence component.
    :param TemplatedSequenceComponent component: a templated sequence component
    :return: component nomenclature representation
    :raises ValueError: if location isn't a SequenceLocation or if unable
    to retrieve region or location
    """
    if component.region and component.region.location:
        location = component.region.location
        if isinstance(location, SequenceLocation):
            sequence_id = location.sequence_id
            start = location.interval.start.number
            end = location.interval.end.number
            return f"{sequence_id}: g.{start}_{end}({component.strand.value})"
        else:
            raise ValueError
    else:
        raise ValueError


def gene_nomenclature(component: GeneComponent) -> str:
    """Return fusion nomenclature for gene component.
    :param GeneComponent component: a gene component object
    :return: component nomenclature representation
    :raises ValueError: if unable to retrieve gene ID
    """
    if component.gene_descriptor.gene_id:
        gene_id = gene_id = component.gene_descriptor.gene_id

    if component.gene_descriptor.gene_id:
        gene_id = component.gene_descriptor.gene_id
    elif component.gene_descriptor.gene \
            and component.gene_descriptor.gene.gene_id:
        gene_id = component.gene_descriptor.gene.gene_id
    else:
        raise ValueError
    return f"{component.gene_descriptor.label}({gene_id})"
