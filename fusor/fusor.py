"""Module for creating and modifying fusion objects."""
from typing import Optional
from biocommons.seqrepo import SeqRepo
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
from ga4gh.vrsatile.pydantic.vrs_model import CURIE
from fusor import SEQREPO_DATA_PATH
from fusor.models import Fusion, GenomicRegionComponent


EXTRA_DATA = [
    "sequence_id",
    "location_id"
]


class FUSOR:
    """Class for creating and modifying fusion objects."""

    def __init__(self, seqrepo_data_path: str = SEQREPO_DATA_PATH):
        """Initialize FUSOR class."""
        self.seqrepo = SeqRepo(seqrepo_data_path)

    def set_fusion_extras(self, fusion: Fusion,
                          set_all: bool = True,
                          fields: Optional[list[str]] = None,
                          target_namespace: str = "ga4gh") -> Fusion:
        """Set extra data to Fusion object.

        :param Fusion fusion: A valid Fusion object
        :param bool set_all: `True` if all extra data fields (see `EXTRA_DATA`)
            will be set in fusion object. `False` if only select fields will
            be provided.
        :param list fields: Select fields that will be set. Must be a subset of
            `EXTRA_DATA`
        :param str target_namespace: The namespace of identifiers to return
            for `sequence_id`. Default is `ga4gh`
        """
        if set_all:
            self.set_sequence_id(fusion, target_namespace)
        else:
            for field in fields:
                if field == "sequence_id":
                    self.set_sequence_id(
                        fusion, target_namespace=target_namespace)
                elif field == "location_id":
                    self.set_location_id(fusion)

        return fusion

    def set_location_id(self, fusion: Fusion) -> None:
        """Set location_id in fusion object.

        :param Fusion fusion: A valid Fusion object
        """
        for transcript_component in fusion.transcript_components:
            if isinstance(transcript_component, GenomicRegionComponent):
                location = transcript_component.region.location
                location_id = \
                    ga4gh_identify(models.Location(**location.dict()))
                transcript_component.region.location_id = location_id

    def set_sequence_id(self, fusion: Fusion,
                        target_namespace: str = "ga4gh") -> None:
        """Set sequence_id in fusion object.

        :param Fusion fusion: A valid Fusion object
        :param str target_namespace: The namespace of identifiers to return
            for `sequence_id`. Default is `ga4gh`
        """
        for transcript_component in fusion.transcript_components:
            if isinstance(transcript_component, GenomicRegionComponent):
                location = transcript_component.region.location
                if location.type == "SequenceLocation":
                    transcript_component.region.location.sequence_id = \
                        self.translate_identifier(location.sequence_id, target_namespace)  # noqa: E501

    def translate_identifier(self, ac: str,
                             target_namespace: str = "ga4gh") -> Optional[CURIE]:
        """Return `target_namespace` identifier for accession provided.

        :param str ac: Identifier accession
        :param str target_namespace: The namespace of identifiers to return.
            Default is `ga4gh`
        :return: Identifier for `target_namespace`
        """
        ga4gh_identifiers = self.seqrepo.translate_identifier(
            ac, target_namespaces=target_namespace)
        if ga4gh_identifiers:
            return ga4gh_identifiers[0]
        return None
