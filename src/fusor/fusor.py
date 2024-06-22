"""Module for modifying fusion objects."""
import logging
import re
from urllib.parse import quote

from biocommons.seqrepo import SeqRepo
from bioutils.accessions import coerce_namespace
from cool_seq_tool.routers import CoolSeqTool
from cool_seq_tool.schemas import ResidueMode
from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models
from ga4gh.vrsatile.pydantic.vrs_models import (
    CURIE,
    Number,
    SequenceInterval,
    SequenceLocation,
    VRSTypes,
)
from ga4gh.vrsatile.pydantic.vrsatile_models import GeneDescriptor, LocationDescriptor
from gene.database import AbstractDatabase as GeneDatabase
from gene.database import create_db
from gene.query import QueryHandler
from pydantic import ValidationError

from fusor.exceptions import FUSORParametersException, IDTranslationException
from fusor.models import (
    AdditionalFields,
    Assay,
    AssayedFusion,
    AssayedFusionElements,
    BaseStructuralElement,
    CategoricalFusion,
    CategoricalFusionElements,
    CausativeEvent,
    DomainStatus,
    Evidence,
    FunctionalDomain,
    Fusion,
    FusionType,
    GeneElement,
    LinkerElement,
    MultiplePossibleGenesElement,
    RegulatoryClass,
    RegulatoryElement,
    Strand,
    StructuralElementType,
    TemplatedSequenceElement,
    TranscriptSegmentElement,
    UnknownGeneElement,
)
from fusor.nomenclature import (
    gene_nomenclature,
    reg_element_nomenclature,
    templated_seq_nomenclature,
    tx_segment_nomenclature,
)
from fusor.tools import translate_identifier

_logger = logging.getLogger(__name__)


class FUSOR:
    """Class for modifying fusion objects."""

    def __init__(
        self,
        seqrepo_data_path: str | None = None,
        gene_database: GeneDatabase | None = None,
        uta_db_url: str | None = None,
    ) -> None:
        """Initialize FUSOR class.

        :param seqrepo_data_path: Path to SeqRepo data directory
        :param gene_database: gene normalizer database instance
        :param uta_db_url: Postgres URL for UTA
        """
        if not gene_database:
            gene_database = create_db()
        self.gene_normalizer = QueryHandler(gene_database)

        cst_params = {}
        if uta_db_url:
            cst_params["db_url"] = uta_db_url
        if seqrepo_data_path:
            cst_params["sr"] = SeqRepo(seqrepo_data_path)
        self.cool_seq_tool = CoolSeqTool(**cst_params)
        self.seqrepo = self.cool_seq_tool.seqrepo_access.sr

    @staticmethod
    def _contains_element_type(kwargs: dict, elm_type: StructuralElementType) -> bool:
        """Check if fusion contains element of a specific type. Helper method
        for inferring fusion type.
        :param Dict kwargs: keyword args given to fusion method
        :param StructuralElementType elm_type: element type to match
        :return: True if at least one element of given type is found,
        False otherwise.
        """
        for c in kwargs["structural_elements"]:
            if (isinstance(c, dict) and c.get("type") == elm_type) or (
                isinstance(c, BaseStructuralElement) and c.type == elm_type
            ):
                return True
        return False

    def fusion(self, fusion_type: FusionType | None = None, **kwargs) -> Fusion:
        """Construct fusion object.

        :param Optional[FusionType] fusion_type: explicitly specify fusion type.
            Unecessary if providing fusion object in keyword args that includes `type`
            attribute.
        :return: constructed fusion object if successful
        :raise: FUSORParametersException if fusion type unable to be determined,
            or if incorrect fusion parameters are provided
        """
        # try explicit type param
        explicit_type = kwargs.get("type")
        if not fusion_type and explicit_type:
            if explicit_type in FusionType.values():
                fusion_type = explicit_type
                kwargs.pop("type")
            else:
                msg = f"Invalid type parameter: {explicit_type}"
                raise FUSORParametersException(msg)
        fusion_fn = None
        if fusion_type:
            if fusion_type == FusionType.CATEGORICAL_FUSION:
                fusion_fn = self.categorical_fusion
            elif fusion_type == FusionType.ASSAYED_FUSION:
                fusion_fn = self.assayed_fusion
            else:
                msg = f"Invalid fusion_type parameter: {fusion_type}"
                raise FUSORParametersException(msg)
        else:
            # try to infer from provided attributes
            categorical_attributes = any(
                [
                    "critical_functional_domains" in kwargs,
                    "r_frame_preserved" in kwargs,
                    self._contains_element_type(
                        kwargs, StructuralElementType.MULTIPLE_POSSIBLE_GENES_ELEMENT
                    ),
                ]
            )
            assayed_attributes = any(
                [
                    "causative_event" in kwargs,
                    "assay" in kwargs,
                    self._contains_element_type(
                        kwargs, StructuralElementType.UNKNOWN_GENE_ELEMENT
                    ),
                ]
            )
            if categorical_attributes and assayed_attributes:
                msg = "Received conflicting attributes"
                raise FUSORParametersException(msg)
            if categorical_attributes and not assayed_attributes:
                fusion_fn = self.categorical_fusion
            elif assayed_attributes and not categorical_attributes:
                fusion_fn = self.assayed_fusion
        if fusion_fn is None:
            msg = "Unable to determine fusion type"
            raise FUSORParametersException(msg)
        try:
            return fusion_fn(**kwargs)
        except TypeError as e:
            msg = f"Unable to construct fusion with provided args: {e}"
            raise FUSORParametersException(msg) from e

    @staticmethod
    def categorical_fusion(
        structural_elements: CategoricalFusionElements,
        regulatory_element: RegulatoryElement | None = None,
        critical_functional_domains: list[FunctionalDomain] | None = None,
        r_frame_preserved: bool | None = None,
    ) -> CategoricalFusion:
        """Construct a categorical fusion object
        :param CategoricalFusionElements structural_elements: elements
            constituting the fusion
        :param Optional[RegulatoryElement] regulatory_element: affected
            regulatory element
        :param Optional[List[FunctionalDomain]] critical_functional_domains: lost or preserved
            functional domains
        :param Optional[bool] r_frame_preserved: `True` if reading frame is
            preserved.  `False` otherwise
        :return: CategoricalFusion if construction successful
        :raise: FUSORParametersException if given incorrect fusion properties
        """
        try:
            fusion = CategoricalFusion(
                structural_elements=structural_elements,
                critical_functional_domains=critical_functional_domains,
                r_frame_preserved=r_frame_preserved,
                regulatory_element=regulatory_element,
            )
        except ValidationError as e:
            raise FUSORParametersException(str(e)) from e
        return fusion

    @staticmethod
    def assayed_fusion(
        structural_elements: AssayedFusionElements,
        causative_event: CausativeEvent | None = None,
        assay: Assay | None = None,
        regulatory_element: RegulatoryElement | None = None,
    ) -> AssayedFusion:
        """Construct an assayed fusion object
        :param AssayedFusionElements structural_elements: elements constituting the
            fusion
        :param Optional[Event] causative_event: event causing the fusion
        :param Optional[Assay] assay: how knowledge of the fusion was obtained
        :param Optional[RegulatoryElement] regulatory_element: affected regulatory
            elements
        :return: Tuple containing optional AssayedFusion if construction successful,
            and any relevant validation warnings
        """
        try:
            fusion = AssayedFusion(
                structural_elements=structural_elements,
                regulatory_element=regulatory_element,
                causative_event=causative_event,
                assay=assay,
            )
        except ValidationError as e:
            raise FUSORParametersException(str(e)) from e
        return fusion

    async def transcript_segment_element(
        self,
        tx_to_genomic_coords: bool = True,
        use_minimal_gene_descr: bool = True,
        seq_id_target_namespace: str | None = None,
        **kwargs,
    ) -> tuple[TranscriptSegmentElement | None, list[str] | None]:
        """Create transcript segment element

        :param bool tx_to_genomic_coords: `True` if going from transcript
            to genomic coordinates. `False` if going from genomic to
            transcript exon coordinates.
        :param bool use_minimal_gene_descr: `True` if minimal gene descriptor
            (`id`, `gene_id`, `label`) will be used. `False` if
            gene-normalizer's gene descriptor will be used
        :param Optional[str] seq_id_target_namespace: If want to use digest for
            `sequence_id`, set this to the namespace you want the digest for.
            Otherwise, leave as `None`.
        :param kwargs:
            If `tx_to_genomic_coords`, possible key word arguments:
                (From cool_seq_tool.transcript_to_genomic_coords)
                gene: Optional[str] = None, transcript: str = None,
                exon_start: Optional[int] = None,
                exon_start_offset: Optional[int] = 0,
                exon_end: Optional[int] = None,
                exon_end_offset: Optional[int] = 0
            else:
                (From cool_seq_tool.genomic_to_transcript_exon_coordinates)
                chromosome: Union[str, int], start: Optional[int] = None,
                end: Optional[int] = None, strand: Optional[int] = None,
                transcript: Optional[str] = None, gene: Optional[str] = None,
                residue_mode: ResidueMode = ResidueMode.RESIDUE
        :return: Transcript Segment Element, warning
        """
        if tx_to_genomic_coords:
            data = await self.cool_seq_tool.ex_g_coords_mapper.transcript_to_genomic_coordinates(
                **kwargs
            )
        else:
            if "chromosome" in kwargs and kwargs.get("chromosome") is None:
                msg = (
                    "`chromosome` is required when going from genomic to"
                    " transcript exon coordinates"
                )
                _logger.warning(msg)
                return None, [msg]
            residue_mode = kwargs.get("residue_mode")
            # TODO: Remove once fixed in cool_seq_tool
            if residue_mode != ResidueMode.INTER_RESIDUE:
                start = kwargs.get("start")
                kwargs["start"] = start - 1 if start is not None else None
                kwargs["residue_mode"] = "inter-residue"
            data = await self.cool_seq_tool.ex_g_coords_mapper.genomic_to_transcript_exon_coordinates(
                **kwargs
            )

        if data.genomic_data is None:
            return None, data.warnings

        genomic_data = data.genomic_data
        genomic_data.transcript = coerce_namespace(genomic_data.transcript)

        normalized_gene_response = self._normalized_gene_descriptor(
            genomic_data.gene, use_minimal_gene_descr=use_minimal_gene_descr
        )
        if not normalized_gene_response[0] and normalized_gene_response[1]:
            return None, [normalized_gene_response[1]]

        return (
            TranscriptSegmentElement(
                transcript=genomic_data.transcript,
                exon_start=genomic_data.exon_start,
                exon_start_offset=genomic_data.exon_start_offset,
                exon_end=genomic_data.exon_end,
                exon_end_offset=genomic_data.exon_end_offset,
                gene_descriptor=normalized_gene_response[0],
                element_genomic_start=self._location_descriptor(
                    genomic_data.start,
                    genomic_data.start + 1,
                    genomic_data.chr,
                    label=genomic_data.chr,
                    seq_id_target_namespace=seq_id_target_namespace,
                )
                if genomic_data.start
                else None,
                element_genomic_end=self._location_descriptor(
                    genomic_data.end,
                    genomic_data.end + 1,
                    genomic_data.chr,
                    label=genomic_data.chr,
                    seq_id_target_namespace=seq_id_target_namespace,
                )
                if genomic_data.end
                else None,
            ),
            None,
        )

    def gene_element(
        self, gene: str, use_minimal_gene_descr: bool = True
    ) -> tuple[GeneElement | None, str | None]:
        """Create gene element

        :param str gene: Gene
        :param bool use_minimal_gene_descr: `True` if minimal gene descriptor
            (`id`, `gene_id`, `label`) will be used. `False` if
            gene-normalizer's gene descriptor will be used
        :return: GeneElement, warning
        """
        gene_descr, warning = self._normalized_gene_descriptor(
            gene, use_minimal_gene_descr=use_minimal_gene_descr
        )
        if not gene_descr:
            return None, warning
        return GeneElement(gene_descriptor=gene_descr), None

    def templated_sequence_element(
        self,
        start: int,
        end: int,
        sequence_id: str,
        strand: Strand,
        label: str | None = None,
        add_location_id: bool = False,
        residue_mode: ResidueMode = ResidueMode.RESIDUE,
        seq_id_target_namespace: str | None = None,
    ) -> TemplatedSequenceElement:
        """Create templated sequence element

        :param int start: Genomic start
        :param int end: Genomic end
        :param str sequence_id: Chromosome accession for sequence
        :param Strand strand: Strand
        :param Optional[str] label: Label for genomic location
        :param bool add_location_id: `True` if `location_id` will be added
            to `region`. `False` otherwise.
        :param ResidueMode residue_mode: Determines coordinate base used.
            Must be one of `residue` or `inter-residue`.
        :param Optional[str] seq_id_target_namespace: If want to use digest for
            `sequence_id`, set this to the namespace you want the digest for.
            Otherwise, leave as `None`.
        :return: Templated Sequence Element
        """
        if residue_mode == ResidueMode.RESIDUE:
            start -= 1

        region = self._location_descriptor(
            start,
            end,
            sequence_id,
            label=label,
            seq_id_target_namespace=seq_id_target_namespace,
        )

        if add_location_id:
            location_id = self._location_id(region.location.model_dump())
            region.location_id = location_id

        return TemplatedSequenceElement(region=region, strand=strand)

    @staticmethod
    def linker_element(
        sequence: str,
        residue_type: CURIE = "SO:0000348",
    ) -> tuple[LinkerElement | None, str | None]:
        """Create linker element

        :param str sequence: Sequence
        :param CURIE residue_type: Sequence Ontology code for residue type of
            `sequence`
        :return: Tuple containing a complete Linker element and None if
            successful, or a None value and warning message if unsuccessful
        """
        try:
            seq = sequence.upper()
            params = {
                "linker_sequence": {
                    "id": f"fusor.sequence:{seq}",
                    "sequence": seq,
                    "residue_type": residue_type,
                }
            }
            return LinkerElement(**params), None
        except ValidationError as e:
            msg = str(e)
            _logger.warning(msg)
            return None, msg

    @staticmethod
    def multiple_possible_genes_element() -> MultiplePossibleGenesElement:
        """Create a MultiplePossibleGenesElement.

        :return: MultiplePossibleGenesElement
        """
        return MultiplePossibleGenesElement()

    @staticmethod
    def unknown_gene_element() -> UnknownGeneElement:
        """Create unknown gene element

        :return: Unknown Gene element
        """
        return UnknownGeneElement()

    def functional_domain(
        self,
        status: DomainStatus,
        name: str,
        functional_domain_id: CURIE,
        gene: str,
        sequence_id: str,
        start: int,
        end: int,
        use_minimal_gene_descr: bool = True,
        seq_id_target_namespace: str | None = None,
    ) -> tuple[FunctionalDomain | None, str | None]:
        """Build functional domain instance.

        :param DomainStatus status: Status for domain.  Must be either `lost`
            or `preserved`
        :param str name: Domain name
        :param CURIE functional_domain_id: Domain ID
        :param str gene: Gene
        :param str sequence_id: protein sequence on which provided coordinates
            are located
        :param int start: start position on sequence
        :param in end: end position on sequence
        :param bool use_minimal_gene_descr: `True` if minimal gene descriptor
            (`id`, `gene_id`, `label`) will be used. `False` if
            gene-normalizer's gene descriptor will be used
        :param Optional[str] seq_id_target_namespace: If want to use digest for
            `sequence_id`, set this to the namespace you want the digest for.
            Otherwise, leave as `None`.
        :return: Tuple with FunctionalDomain and None value for warnings if
            successful, or a None value and warning message if unsuccessful
        """
        sequence_id_lower = sequence_id.lower()
        if not (sequence_id_lower.startswith("np_")) or (
            sequence_id_lower.startswith("ensp")
        ):
            msg = "Sequence_id must be a protein accession."
            _logger.warning(msg)
            return None, msg

        seq, warning = self.cool_seq_tool.seqrepo_access.get_reference_sequence(
            sequence_id, start, end
        )

        if not seq:
            return None, warning

        gene_descr, warning = self._normalized_gene_descriptor(
            gene, use_minimal_gene_descr=use_minimal_gene_descr
        )
        if not gene_descr:
            return None, warning

        loc_descr = self._location_descriptor(
            start, end, sequence_id, seq_id_target_namespace=seq_id_target_namespace
        )

        try:
            return (
                FunctionalDomain(
                    _id=functional_domain_id,
                    label=name,
                    status=status,
                    associated_gene=gene_descr,
                    sequence_location=loc_descr,
                ),
                None,
            )
        except ValidationError as e:
            msg = str(e)
            _logger.warning(msg)
            return None, msg

    def regulatory_element(
        self,
        regulatory_class: RegulatoryClass,
        gene: str,
        use_minimal_gene_descr: bool = True,
    ) -> tuple[RegulatoryElement | None, str | None]:
        """Create RegulatoryElement
        :param RegulatoryClass regulatory_class: one of {"promoter", "enhancer"}
        :param str gene: gene term to fetch normalized descriptor for
        :param bool use_minimal_gene_descr: whether to use the minimal gene descriptor
        :return: Tuple with RegulatoryElement instance and None value for warnings if
            successful, or a None value and warning message if unsuccessful
        """
        gene_descr, warning = self._normalized_gene_descriptor(
            gene, use_minimal_gene_descr=use_minimal_gene_descr
        )
        if not gene_descr:
            return None, warning

        try:
            return (
                RegulatoryElement(
                    regulatory_class=regulatory_class, associated_gene=gene_descr
                ),
                None,
            )
        except ValidationError as e:
            msg = str(e)
            _logger.warning(msg)
            return None, msg

    def _location_descriptor(
        self,
        start: int,
        end: int,
        sequence_id: str,
        label: str | None = None,
        seq_id_target_namespace: str | None = None,
        use_location_id: bool = False,
    ) -> LocationDescriptor:
        """Create location descriptor

        :param int start: Start position
        :param int end: End position
        :param str sequence_id: Accession for sequence
        :param str label: label for location. If `None`, `sequence_id` will be used as
            Location Descriptor's `id` Else, label will be used as Location
            Descriptor's `id`.
        :param str seq_id_target_namespace: If want to use digest for `sequence_id`,
            set this to the namespace you want the digest for. Otherwise, leave as
            `None`.
        :param bool use_location_id: Takes precedence over `label` or `sequence_id`
            becoming Location Descriptor's id. `True` if  use ga4gh digest as Location
            Descriptor's id. `False`, use default of `label` > `sequence_id`
        """
        seq_id_input = sequence_id

        try:
            sequence_id = coerce_namespace(sequence_id)
        except ValueError:
            if not re.match(CURIE.__metadata__[0].pattern, sequence_id):
                sequence_id = f"sequence.id:{sequence_id}"

        if seq_id_target_namespace:
            try:
                seq_id = translate_identifier(
                    self.seqrepo, sequence_id, target_namespace=seq_id_target_namespace
                )
            except IDTranslationException:
                _logger.warning(
                    "Unable to translate %s using %s as the target namespace",
                    sequence_id,
                    seq_id_target_namespace,
                )
            else:
                sequence_id = seq_id

        location = SequenceLocation(
            sequence_id=sequence_id,
            interval=SequenceInterval(start=Number(value=start), end=Number(value=end)),
        )

        if use_location_id:
            _id = self._location_id(location.model_dump())
        else:
            quote_id = quote(label) if label else quote(seq_id_input)
            _id = f"fusor.location_descriptor:{quote_id}"

        location_descr = LocationDescriptor(id=_id, location=location)

        if label:
            location_descr.label = label
        return location_descr

    def add_additional_fields(
        self,
        fusion: Fusion,
        add_all: bool = True,
        fields: list[AdditionalFields] | None = None,
        target_namespace: str = "ga4gh",
    ) -> Fusion:
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
            self.add_translated_sequence_id(fusion, target_namespace)
            self.add_location_id(fusion)
        else:
            if fields:
                for field in fields:
                    if field == AdditionalFields.SEQUENCE_ID.value:
                        self.add_translated_sequence_id(
                            fusion, target_namespace=target_namespace
                        )
                    elif field == AdditionalFields.LOCATION_ID.value:
                        self.add_location_id(fusion)
                    else:
                        _logger.warning("Invalid field: %s", field)
        return fusion

    def add_location_id(self, fusion: Fusion) -> Fusion:
        """Add `location_id` in fusion object.

        :param Fusion fusion: A valid Fusion object.
        :return: Updated fusion with `location_id` fields set
        """
        for structural_element in fusion.structural_elements:
            if isinstance(structural_element, TemplatedSequenceElement):
                location = structural_element.region.location
                location_id = self._location_id(location.model_dump())
                structural_element.region.location_id = location_id
            elif isinstance(structural_element, TranscriptSegmentElement):
                for element_genomic in [
                    structural_element.element_genomic_start,
                    structural_element.element_genomic_end,
                ]:
                    if element_genomic:
                        location = element_genomic.location
                        if location.type == VRSTypes.SEQUENCE_LOCATION.value:
                            location_id = self._location_id(location.model_dump())
                            element_genomic.location_id = location_id
        if isinstance(fusion, CategoricalFusion) and fusion.critical_functional_domains:
            for domain in fusion.critical_functional_domains:
                location = domain.sequence_location.location
                location_id = self._location_id(location.model_dump())
                domain.sequence_location.location_id = location_id
        if fusion.regulatory_element:
            element = fusion.regulatory_element
            if element.feature_location:
                location = element.feature_location
                if location.type == VRSTypes.SEQUENCE_LOCATION.value:
                    location_id = self._location_id(location.model_dump())
                    element.feature_location.location_id = location_id
        return fusion

    @staticmethod
    def _location_id(location: dict) -> CURIE:
        """Return GA4GH digest for location

        :param dict location: VRS Location represented as a dict
        :return: GA4GH digest
        """
        return ga4gh_identify(models.Location(**location))

    def add_translated_sequence_id(
        self, fusion: Fusion, target_namespace: str = "ga4gh"
    ) -> Fusion:
        """Translate sequence_ids in fusion object.

        :param Fusion fusion: A valid Fusion object
        :param str target_namespace: ID namespace to translate sequence IDs to
        :return: Updated fusion with `sequence_id` fields set
        """
        for element in fusion.structural_elements:
            if isinstance(element, TemplatedSequenceElement):
                location = element.region.location
                if location.type == VRSTypes.SEQUENCE_LOCATION.value:
                    try:
                        new_id = translate_identifier(
                            self.seqrepo, location.sequence_id, target_namespace
                        )
                    except IDTranslationException:
                        pass
                    else:
                        element.region.location.sequence_id = new_id
            elif isinstance(element, TranscriptSegmentElement):
                for loc_descr in [
                    element.element_genomic_start,
                    element.element_genomic_end,
                ]:
                    if loc_descr:
                        location = loc_descr.location
                        if location.type == VRSTypes.SEQUENCE_LOCATION.value:
                            try:
                                new_id = translate_identifier(
                                    self.seqrepo, location.sequence_id, target_namespace
                                )
                            except IDTranslationException:
                                continue
                            loc_descr.location.sequence_id = new_id
        if fusion.type == "CategoricalFusion" and fusion.critical_functional_domains:
            for domain in fusion.critical_functional_domains:
                if (
                    domain.sequence_location
                    and domain.sequence_location.location
                    and (domain.sequence_location.location.type == "SequenceLocation")
                ):
                    try:
                        new_id = translate_identifier(
                            self.seqrepo,
                            domain.sequence_location.location.sequence_id,
                            target_namespace,
                        )
                    except IDTranslationException:
                        continue
                    domain.sequence_location.location.sequence_id = new_id
        return fusion

    def add_gene_descriptor(self, fusion: Fusion) -> Fusion:
        """Add additional fields to `gene_descriptor` in fusion object

        :param Fusion fusion: A valid Fusion object
        :return: Updated fusion with additional fields set in `gene_descriptor`
        """
        properties = [fusion.structural_elements]
        if fusion.type == FusionType.CATEGORICAL_FUSION:
            properties.append(fusion.critical_functional_domains)

        for prop in properties:
            for obj in prop:
                if "gene_descriptor" in obj.model_fields:
                    label = obj.gene_descriptor.label
                    norm_gene_descr, _ = self._normalized_gene_descriptor(
                        label, use_minimal_gene_descr=False
                    )
                    if norm_gene_descr:
                        obj.gene_descriptor = norm_gene_descr
        if fusion.regulatory_element and fusion.regulatory_element.associated_gene:
            reg_el = fusion.regulatory_element
            label = reg_el.associated_gene.label
            norm_gene_descr, _ = self._normalized_gene_descriptor(
                label, use_minimal_gene_descr=False
            )
            if norm_gene_descr:
                reg_el.associated_gene = norm_gene_descr
        return fusion

    def _normalized_gene_descriptor(
        self, query: str, use_minimal_gene_descr: bool = True
    ) -> tuple[GeneDescriptor | None, str | None]:
        """Return gene descriptor from normalized response.

        :param str query: Gene query
        :param bool use_minimal_gene_descr: `True` if minimal gene descriptor
            (`id`, `gene_id`, `label`) will be used. `False` if
            gene-normalizer's gene descriptor will be used
        :return: Tuple with gene descriptor and None value for warnings if
            successful, and None value with warning string if unsuccessful
        """
        gene_norm_resp = self.gene_normalizer.normalize(query)
        if gene_norm_resp.match_type:
            gene_descr = gene_norm_resp.gene_descriptor
            if use_minimal_gene_descr:
                gene_descr = GeneDescriptor(
                    id=gene_descr.id, gene_id=gene_descr.gene_id, label=gene_descr.label
                )
            return gene_descr, None
        return None, f"gene-normalizer unable to normalize {query}"

    def generate_nomenclature(self, fusion: Fusion) -> str:
        """Generate human-readable nomenclature describing provided fusion
        :param Fusion fusion: a valid fusion
        :return: string summarizing fusion in human-readable way per
            VICC fusion curation nomenclature
        """
        parts = []
        element_genes = []
        if fusion.regulatory_element:
            parts.append(
                reg_element_nomenclature(fusion.regulatory_element, self.seqrepo)
            )
        for element in fusion.structural_elements:
            if isinstance(element, MultiplePossibleGenesElement):
                parts.append("v")
            elif isinstance(element, UnknownGeneElement):
                parts.append("?")
            elif isinstance(element, LinkerElement):
                parts.append(element.linker_sequence.sequence)
            elif isinstance(element, TranscriptSegmentElement):
                if not any(
                    [gene == element.gene_descriptor.label for gene in element_genes]  # noqa: C419
                ):
                    parts.append(tx_segment_nomenclature(element))
            elif isinstance(element, TemplatedSequenceElement):
                parts.append(templated_seq_nomenclature(element, self.seqrepo))
            elif isinstance(element, GeneElement):
                if not any(
                    [gene == element.gene_descriptor.label for gene in element_genes]  # noqa: C419
                ):
                    parts.append(gene_nomenclature(element))
            else:
                raise ValueError
        if (
            isinstance(fusion, AssayedFusion)
            and fusion.assay
            and fusion.assay.fusion_detection == Evidence.INFERRED
        ):
            divider = "(::)"
        else:
            divider = "::"
        return divider.join(parts)
