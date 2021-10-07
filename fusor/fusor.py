"""Module for modifying fusion objects."""
from typing import Optional, List
from biocommons.seqrepo import SeqRepo
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
from ga4gh.vrsatile.pydantic.vrs_model import CURIE, VRSTypes
from fusor import SEQREPO_DATA_PATH
from fusor.models import Fusion, TemplatedSequenceComponent,\
    AdditionalFields, TranscriptSegmentComponent
from fusor import logger


class FUSOR:
    """Class for modifying fusion objects."""

    def __init__(self, seqrepo_data_path: str = SEQREPO_DATA_PATH) -> None:
        """Initialize FUSOR class.

        :param str seqrepo_data_path: Path to SeqRepo data directory
        """
        self.seqrepo = SeqRepo(seqrepo_data_path)

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

    def add_location_id(self, fusion: Fusion) -> None:
        """Add `location_id` in fusion object.

        :param Fusion fusion: A valid Fusion object
        """
        for structural_component in fusion.structural_components:
            if isinstance(structural_component, TemplatedSequenceComponent):
                location = structural_component.region.location
                location_id = \
                    ga4gh_identify(models.Location(**location.dict()))
                structural_component.region.location_id = location_id
            elif isinstance(structural_component, TranscriptSegmentComponent):
                location = structural_component.component_genomic_region.location  # noqa: E501
                location_id = \
                    ga4gh_identify(models.Location(**location.dict()))
                if location.type == VRSTypes.SEQUENCE_LOCATION.value:
                    structural_component.component_genomic_region.location_id = location_id  # noqa: E501

    def add_sequence_id(self, fusion: Fusion,
                        target_namespace: str = "ga4gh") -> None:
        """Add sequence_id in fusion object.

        :param Fusion fusion: A valid Fusion object
        :param str target_namespace: The namespace of identifiers to return
            for `sequence_id`. Default is `ga4gh`
        """
        for structural_component in fusion.structural_components:
            if isinstance(structural_component, TemplatedSequenceComponent):
                location = structural_component.region.location
                if location.type == VRSTypes.SEQUENCE_LOCATION.value:
                    structural_component.region.location.sequence_id = \
                        self.translate_identifier(location.sequence_id, target_namespace)  # noqa: E501
            elif isinstance(structural_component, TranscriptSegmentComponent):
                location = structural_component.component_genomic_region.location  # noqa: #501
                if location.type == VRSTypes.SEQUENCE_LOCATION.value:
                    structural_component.component_genomic_region.\
                        location.sequence_id = self.translate_identifier(location.sequence_id, target_namespace)  # noqa: E501

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
