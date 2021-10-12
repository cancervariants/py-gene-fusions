"""Module for modifying fusion objects."""
from typing import Optional, List, Union, Tuple, Dict
from biocommons.seqrepo import SeqRepo
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
from ga4gh.vrsatile.pydantic.vrs_model import CURIE, VRSTypes, \
    SequenceLocation, Number, SequenceInterval
from ga4gh.vrsatile.pydantic.vrsatile_model import GeneDescriptor, \
    SequenceDescriptor, LocationDescriptor
from pydantic.error_wrappers import ValidationError
from uta_tools.uta_tools import UTATools

from fusor import SEQREPO_DATA_PATH
from gene.query import QueryHandler
from fusor.models import Fusion, TemplatedSequenceComponent, \
    AdditionalFields, TranscriptSegmentComponent, GeneComponent, \
    LinkerComponent, UnknownGeneComponent, RegulatoryElement, Event, \
    DomainStatus, CriticalDomain, Strand
from fusor import logger, UTA_DB_URL
from bioutils.accessions import coerce_namespace
from uta_tools.schemas import ResidueMode
from urllib.parse import quote


class FUSOR:
    """Class for modifying fusion objects."""

    def __init__(self,
                 seqrepo_data_path: str = SEQREPO_DATA_PATH,
                 dynamodb_url: str = "",
                 dynamodb_region: str = "us-east-2",
                 db_url: str = UTA_DB_URL, db_pwd: str = "",
                 ) -> None:
        """Initialize FUSOR class.

        :param str seqrepo_data_path: Path to SeqRepo data directory
        :param str dynamodb_url: URL to gene-normalizer database source.
            Can also set environment variable `GENE_NORM_DB_URL`.
        :param str dynamodb_region: AWS default region for gene-normalizer.
        """
        self.seqrepo = SeqRepo(seqrepo_data_path)
        self.gene_normalizer = QueryHandler(
            db_url=dynamodb_url, db_region=dynamodb_region)
        self.uta_tools = UTATools(db_url=db_url, db_pwd=db_pwd)

    def fusion(
            self,
            structural_components: List[Union[
                TranscriptSegmentComponent, GeneComponent,
                TemplatedSequenceComponent, LinkerComponent,
                UnknownGeneComponent]],
            r_frame_preserved: Optional[bool] = None,
            causative_event: Optional[Event] = None,
            regulatory_elements: Optional[RegulatoryElement] = None) -> Fusion:
        return Fusion(r_frame_preserved=r_frame_preserved,
                      tructural_components=structural_components,
                      causative_event=causative_event,
                      regulatory_elements=regulatory_elements)

    async def transcript_to_genomic_coordinates(
            self, gene: Optional[str] = None, transcript: str = None,
            exon_start: Optional[int] = None, exon_start_offset: Optional[int] = 0,  # noqa: E501
            exon_end: Optional[int] = None, exon_end_offset: Optional[int] = 0,
            **kwargs):
        pass

    async def genomic_to_transcript_exon_coordinates(
            self, chromosome: Union[str, int], start: Optional[int] = None,
            end: Optional[int] = None, strand: Optional[int] = None,
            transcript: Optional[str] = None, gene: Optional[str] = None,
            residue_mode: ResidueMode = ResidueMode.RESIDUE,
            **kwargs):
        pass

    async def transcript_segment_component(
            self, use_exon_coords: bool = True,
            use_minimal_gene_descr: bool = True, **kwargs
    ) -> Tuple[Optional[TranscriptSegmentComponent], Optional[str]]:
        if use_exon_coords:
            data = await self.uta_tools.transcript_to_genomic_coordinates(**kwargs)
        else:
            if "chromosome" in kwargs and kwargs.get("chromosome") is None:
                msg = "`chromosome` is required when going from genomic to" \
                      " transcript exon coordinates"
                logger.warning(msg)
                return None, msg
            data = await self.uta_tools.genomic_to_transcript_exon_coordinates(**kwargs)

        if data.genomic_data is None:
            return None, data.warnings

        genomic_data = data.genomic_data
        genomic_data.transcript = coerce_namespace(genomic_data.transcript)

        return TranscriptSegmentComponent(
            transcript=genomic_data.transcript,
            exon_start=genomic_data.exon_start,
            exon_start_offset=genomic_data.exon_start_offset,
            exon_end=genomic_data.exon_end,
            exon_end_offset=genomic_data.exon_end_offset,
            gene_descriptor=self._normalized_gene_descriptor(
                genomic_data.gene,
                use_minimal_gene_descr=use_minimal_gene_descr)[0],
            component_genomic_start=self._location_descriptor(
                genomic_data.start, genomic_data.start + 1, genomic_data.chr,
                label=genomic_data.chr,
                seq_id_target_namespace="ga4gh") if genomic_data.start else None,  # noqa: E501
            component_genomic_end=self._location_descriptor(
                genomic_data.end - 1, genomic_data.end, genomic_data.chr,
                label=genomic_data.chr,
                seq_id_target_namespace="ga4gh") if genomic_data.end else None,
        ), None

    def gene_component(self, gene: str,
                       use_minimal_gene_descr: bool = True) -> Tuple[Optional[GeneComponent], Optional[str]]:
        gene_descr, warning = self._normalized_gene_descriptor(
            gene, use_minimal_gene_descr=use_minimal_gene_descr)
        if not gene_descr:
            return None, warning
        else:
            return gene_descr, None

    def templated_sequence_component(
            self, start: int, end: int, sequence_id: str, strand: Strand,
            label: Optional[str] = None, add_location_id: bool = False,
            residue_mode: ResidueMode = ResidueMode.RESIDUE,
            seq_id_target_namespace: Optional[str] = "ga4gh") -> TemplatedSequenceComponent:

        if residue_mode == ResidueMode.RESIDUE:
            start -= 1

        region = self._location_descriptor(
            start, end, sequence_id, label=label,
            seq_id_target_namespace=seq_id_target_namespace)

        if add_location_id:
            location_id = self._location_id(region.location)
            region.location_id = location_id

        return TemplatedSequenceComponent(region=region, strand=strand)

    def linker_component(
            self, sequence: str,
            residue_type: CURIE = "SO:0000348") -> Tuple[Optional[LinkerComponent], Optional[str]]:
        try:
            sequence_descriptor = SequenceDescriptor(
                id=f"fusor.sequence:{sequence.upper()}",
                sequence=sequence,
                residue_type=residue_type
            )
        except ValidationError as e:
            msg = str(e)
            logger.warning(msg)
            return None, msg
        else:
            return LinkerComponent(linker_sequence=sequence_descriptor), None

    def unknown_gene_component(self) -> UnknownGeneComponent:
        return UnknownGeneComponent()

    def critical_domain(self, status: DomainStatus, name: str,
                        critical_domain_id: CURIE, gene: str,
                        use_minimal_gene_descr: bool = True) -> Tuple[Optional[CriticalDomain], Optional[str]]:
        gene_descr, warning = self._normalized_gene_descriptor(
            gene, use_minimal_gene_descr=use_minimal_gene_descr)
        if not gene_descr:
            return None, warning

        try:
            return CriticalDomain(
                id=critical_domain_id,
                name=name,
                status=status,
                gene_descriptor=gene_descr
            ), None
        except ValidationError as e:
            msg = str(e)
            logger.warning(msg)
            return None, msg

    def _location_descriptor(self, start: int, end: int, sequence_id: str,
                             label: Optional[str] = None,
                             seq_id_target_namespace: Optional[str] = None):
        try:
            sequence_id = coerce_namespace(sequence_id)
        except ValueError:
            try:
                CURIE(__root__=sequence_id)
            except ValidationError:
                sequence_id = f"sequence.id:{sequence_id}"

        if seq_id_target_namespace:
            seq_id = self.translate_identifier(
                sequence_id, target_namespace=seq_id_target_namespace)
            if seq_id:
                sequence_id = seq_id
            else:
                logger.warning(f"Unable to translate {sequence_id} using"
                               f" {seq_id_target_namespace} as the target"
                               f" namespace")

        location = SequenceLocation(
            sequence_id=sequence_id,
            interval=SequenceInterval(start=Number(value=start),
                                      end=Number(value=end))
        )

        location_descr = LocationDescriptor(
            id=f"fusor.location_descriptor:{quote(label)}",
            location=location
        )
        if label:
            location_descr.label = label
        return location_descr

    def add_additional_fields(self, fusion: Fusion,
                              add_all: bool = True,
                              fields: Optional[List[AdditionalFields]] = None,
                              target_namespace: str = "ga4gh") -> Fusion:
        """Add additional fields to Fusion object.
        Possible fields are shown in `AdditionalFields`

        :param Fusion fusion: A valid Fusion object
        :param bool add_all: `True` if all additional fields  will be added
            in fusion object. `False` if only select fields will be provided.
            If set to `True`, will always take precedence over `fields`.
        :param list fields: Select fields that will be set. Must be a subset of
            `AdditionalFields`
        :param str target_namespace: The namespace of identifiers to return
            for `sequence_id`. Default is `ga4gh`
        :return: Updated fusion with specified fields set
        """
        if add_all:
            self.add_sequence_id(fusion, target_namespace)
            self.add_location_id(fusion)
        else:
            for field in fields:
                if field == AdditionalFields.SEQUENCE_ID.value:
                    self.add_sequence_id(
                        fusion, target_namespace=target_namespace)
                elif field == AdditionalFields.LOCATION_ID.value:
                    self.add_location_id(fusion)
                else:
                    logger.warning(f"Invalid field: {field}")

        return fusion

    def add_location_id(self, fusion: Fusion) -> Fusion:
        """Add `location_id` in fusion object.

        :param Fusion fusion: A valid Fusion object
        :return: Updated fusion with `location_id` fields set
        """
        for structural_component in fusion.structural_components:
            if isinstance(structural_component, TemplatedSequenceComponent):
                location = structural_component.region.location
                location_id = self._location_id(location.dict())
                structural_component.region.location_id = location_id
            elif isinstance(structural_component, TranscriptSegmentComponent):
                for component_genomic in [
                    structural_component.component_genomic_start,
                    structural_component.component_genomic_end
                ]:
                    if component_genomic:
                        location = component_genomic.location
                        if location.type == VRSTypes.SEQUENCE_LOCATION.value:
                            location_id = self._location_id(location.dict())
                            component_genomic.location_id = location_id
        return fusion

    def _location_id(self, location: Dict):
        return ga4gh_identify(models.Location(location))

    def add_sequence_id(self, fusion: Fusion,
                        target_namespace: str = "ga4gh") -> Fusion:
        """Add sequence_id in fusion object.

        :param Fusion fusion: A valid Fusion object
        :param str target_namespace: The namespace of identifiers to return
            for `sequence_id`. Default is `ga4gh`
        :return: Updated fusion with `sequence_id` fields set
        """
        for structural_component in fusion.structural_components:
            if isinstance(structural_component, TemplatedSequenceComponent):
                location = structural_component.region.location
                if location.type == VRSTypes.SEQUENCE_LOCATION.value:
                    structural_component.region.location.sequence_id = \
                        self.translate_identifier(location.sequence_id, target_namespace)  # noqa: E501
            elif isinstance(structural_component, TranscriptSegmentComponent):
                for component_genomic in [
                    structural_component.component_genomic_start,
                    structural_component.component_genomic_end
                ]:
                    if component_genomic:
                        location = component_genomic.location
                        if location.type == VRSTypes.SEQUENCE_LOCATION.value:
                            component_genomic.location.sequence_id = \
                                self.translate_identifier(location.sequence_id, target_namespace)  # noqa: E501
        return fusion

    def add_gene_descriptor(self, fusion: Fusion) -> Fusion:
        """Add additional fields to `gene_descriptor` in fusion object

        :param Fusion fusion: A valid Fusion object
        :return: Updated fusion with additional fields set in `gene_descriptor`
        """
        for field in [fusion.protein_domains, fusion.structural_components,
                      fusion.regulatory_elements]:
            for obj in field:
                if "gene_descriptor" in obj.__fields__.keys():
                    norm_gene_descr, _ = \
                        self._normalized_gene_descriptor(obj.gene_descriptor.label, use_minimal_gene_descr=False)  # noqa: E501
                    if norm_gene_descr:
                        obj.gene_descriptor = norm_gene_descr
        return fusion

    def _normalized_gene_descriptor(
            self, query: str,
            use_minimal_gene_descr: bool = True) -> Tuple[Optional[GeneDescriptor], Optional[str]]:
        """Return gene descriptor from normalized response.

        :param str query: Gene query
        :return: Gene Descriptor for query
        """
        gene_norm_resp = self.gene_normalizer.normalize(query)
        if gene_norm_resp.match_type:
            gene_descr = gene_norm_resp.gene_descriptor
            if use_minimal_gene_descr:
                gene_descr = GeneDescriptor(
                    id=gene_descr.id,
                    gene_id=gene_descr.id,
                    label=gene_descr.label
                )
            return gene_descr, None
        else:
            return None, f"gene-normalizer unable to normalize {query}"

    def translate_identifier(self, ac: str,
                             target_namespace: str = "ga4gh") -> Optional[CURIE]:  # noqa: E501
        """Return `target_namespace` identifier for accession provided.

        :param str ac: Identifier accession
        :param str target_namespace: The namespace of identifiers to return.
            Default is `ga4gh`
        :return: Identifier for `target_namespace`
        """
        try:
            ga4gh_identifiers = self.seqrepo.translate_identifier(
                ac, target_namespaces=target_namespace)
        except KeyError as e:
            logger.warning(f"Unable to get translated identifier: {e}")
            return None

        if ga4gh_identifiers:
            return ga4gh_identifiers[0]
        return None
