"""Schemas for fusion callers used in translator.py"""

from enum import Enum
from typing import Literal

from pydantic import BaseModel, Field


class BaseModelForbidExtra(BaseModel, extra="forbid"):
    """Base Pydantic model class with extra values forbidden."""


class FusionCallerTypes(str, Enum):
    """Define FusionCaller type values"""

    JAFFA = "JAFFA"


class JAFFA(BaseModel):
    """Define parameters for JAFFA model"""

    type: Literal[FusionCallerTypes.JAFFA] = FusionCallerTypes.JAFFA
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
