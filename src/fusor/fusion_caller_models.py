"""Schemas for fusion callers used in translator.py"""

from enum import Enum
from typing import Literal

from pydantic import BaseModel, Field


class BaseModelForbidExtra(BaseModel, extra="forbid"):
    """Base Pydantic model class with extra values forbidden."""


class Caller(str, Enum):
    """Define different supported callers"""

    JAFFA = "JAFFA"
    STAR_FUSION = "STAR-Fusion"
    FUSION_CATCHER = "FusionCatcher"
    FUSION_MAP = "FusionMap"
    ARRIBA = "Arriba"
    CICERO = "CICERO"
    MAPSPLICE = "MapSplice"
    ENFUSION = "EnFusion"
    GENIE = "GENIE"


class JAFFA(BaseModel):
    """Define parameters for JAFFA model"""

    type: Literal[Caller.JAFFA] = Caller.JAFFA
    fusion_genes: str = Field(
        ..., description="A string containing the two fusion partners"
    )
    chrom1: str = Field(
        ..., description="The chromosome indicated in the chrom1 column"
    )
    base1: int = Field(
        ..., description="The genomic position indicated in the base1 column"
    )
    chrom2: str = Field(
        ..., description="The chromosome indicated in the chrom2 column"
    )
    base2: int = Field(
        ..., description="The genomic position indicated in the base2 column"
    )
    rearrangement: bool = Field(
        ..., description=" A boolean indicating if a rearrangement occured"
    )
    classification: str = Field(
        ..., description="The classification associated with the called fusion"
    )
    inframe: bool = Field(
        ..., description="A boolean indicating if the fusion occurred in-frame"
    )
    spanning_reads: int = Field(
        ...,
        description="The number of deteced reads that span the junction bewtween the two transcript. Although described as spanning reads, this aligns with our defintion of split reads i.e. reads that have sequence belonging to the fusion partners",
    )
    spanning_pairs: int = Field(
        ...,
        description="The number of detected reads that align entirely on either side of the breakpoint",
    )


class STARFusion(BaseModel):
    """Define parameters for STAR-Fusion model"""

    type: Literal[Caller.STAR_FUSION] = Caller.STAR_FUSION
    left_gene: str = Field(..., description="The gene indicated in the LeftGene column")
    right_gene: str = Field(
        ..., description="The gene indicated in the RightGene column"
    )
    left_breakpoint: str = Field(
        ..., description="The gene indicated in the LeftBreakpoint column"
    )
    right_breakpoint: str = Field(
        ..., description="The gene indicated in the RightBreakpoint column"
    )
    annots: str = Field(..., description="The annotations associated with the fusion")
    junction_read_count: int = Field(
        ...,
        description="The number of RNA-seq fragments that split the junction between the two transcript segments (from STAR-Fusion documentation)",
    )
    spanning_frag_count: int = Field(
        ...,
        description="The number of RNA-seq fragments that encompass the fusion junction such that one read of the pair aligns to a different gene than the other paired-end read of that fragment (from STAR-Fusion documentation)",
    )


class FusionCatcher(BaseModel):
    """Define parameters for FusionCatcher model"""

    type: Literal[Caller.FUSION_CATCHER] = Caller.FUSION_CATCHER
    five_prime_partner: str = Field(
        ..., description="Gene symbol for the 5' fusion partner"
    )
    three_prime_partner: str = Field(
        ..., description="Gene symbol for the 3' fusion partner"
    )
    five_prime_fusion_point: str = Field(
        ...,
        description="Chromosomal position for the 5' end of the fusion junction. This coordinate is 1-based",
    )
    three_prime_fusion_point: str
    predicted_effect: str
    spanning_unique_reads: int
    spanning_reads: int
    fusion_sequence: str
