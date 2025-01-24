"""Schemas for outputs provided by different fusion callers"""

import csv
from abc import ABC, abstractmethod
from enum import Enum
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field


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


class FusionCaller(ABC, BaseModel):
    """ABC for fusion callers"""

    type: Caller
    model_config = ConfigDict(extra="allow")

    @staticmethod
    def _does_file_exist(path: Path) -> None:
        """Check if fusions file exists

        :param path: The path to the file
        :return None
        :raise ValueError if the file does not exist at the specified path
        """
        if not path.exists():
            statement = f"{path!s} does not exist"
            raise ValueError(statement)
        return

    @classmethod
    def _process_fusion_caller_rows(
        cls,
        path: Path,
        column_rename: dict,
        delimeter: str,
    ) -> list["FusionCaller"]:
        """Convert rows of fusion caller output to Pydantic classes

        :param path: The path to the fusions file
        :param column_rename: A dictionary of column mappings
        :param delimeter: The delimeter for the fusions file
        :return: A list of fusions, represented as Pydantic objects
        """
        cls._does_file_exist(path)
        fusions_list = []
        with path.open() as csvfile:
            reader = csv.DictReader(csvfile, delimiter=delimeter)
            for row in reader:
                row = {column_rename.get(key, key): value for key, value in row.items()}
                fusions_list.append(cls(**row))
        return fusions_list

    @abstractmethod
    def load_records(self, path: Path) -> list[Caller]:
        """Abstract method to load records from a fusion caller file."""


class JAFFA(FusionCaller):
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
        ..., description=" A boolean indicating if a rearrangement occurred"
    )
    classification: str = Field(
        ..., description="The classification associated with the called fusion"
    )
    inframe: bool | str = Field(
        ...,
        description="A boolean or string indicating if the fusion occurred in-frame",
    )
    spanning_reads: int = Field(
        ...,
        description="The number of detected reads that span the junction between the two transcript. Although described as spanning reads, this aligns with our definition of split reads i.e. reads that have sequence belonging to the two fusion partners",
    )
    spanning_pairs: int = Field(
        ...,
        description="The number of detected reads that align entirely on either side of the breakpoint",
    )

    @classmethod
    def load_records(cls, path: Path) -> list["JAFFA"]:
        """Load fusions from JAFFA csv file

        :param path: The path to the file of JAFFA fusions
        :return A list of JAFFA objects, or None if the specified file does not exist
        """
        column_rename = {
            "fusion genes": "fusion_genes",
            "spanning reads": "spanning_reads",
            "spanning pairs": "spanning_pairs",
        }
        return cls._process_fusion_caller_rows(path, column_rename, ",")


class STARFusion(FusionCaller):
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

    @classmethod
    def load_records(cls, path: Path) -> list["STARFusion"]:
        """Load fusions from STAR-Fusion tsv file

        :param path: The path to the file of STAR-Fusion fusions
        :return A list of STAR-Fusion objects, or None if the specified file does not exist
        """
        column_rename = {
            "LeftGene": "left_gene",
            "RightGene": "right_gene",
            "LeftBreakpoint": "left_breakpoint",
            "RightBreakpoint": "right_breakpoint",
            "JunctionReadCount": "junction_read_count",
            "SpanningFragCount": "spanning_frag_count",
        }
        return cls._process_fusion_caller_rows(path, column_rename, "\t")


class FusionCatcher(FusionCaller):
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
    three_prime_fusion_point: str = Field(
        ...,
        description="Chromosomal position for the 3' end of the fusion junction. This coordinate is 1-based",
    )
    predicted_effect: str = Field(
        ...,
        description="The predicted effect of the fusion event, created using annotation from the Ensembl database",
    )
    spanning_unique_reads: int = Field(
        ..., description="The number of unique reads that map on the fusion junction"
    )
    spanning_reads: int = Field(
        ..., description="The number of paired reads that support the fusion"
    )
    fusion_sequence: str = Field(
        ..., description="The inferred sequence around the fusion junction"
    )

    @classmethod
    def load_records(cls, path: Path) -> list["FusionCatcher"]:
        """Load fusions from FusionCatcher txt file

        :param path: The path to the file of FusionCatcher fusions
        :return A list of FusionCatcher objects, or None if the specified file does not exist
        """
        column_rename = {
            "Gene_1_symbol(5end_fusion_partner)": "five_prime_partner",
            "Gene_2_symbol(3end_fusion_partner)": "three_prime_partner",
            "Fusion_point_for_gene_1(5end_fusion_partner)": "five_prime_fusion_point",
            "Fusion_point_for_gene_2(3end_fusion_partner)": "three_prime_fusion_point",
            "Predicted_effect": "predicted_effect",
            "Spanning_unique_reads": "spanning_unique_reads",
            "Spanning_pairs": "spanning_reads",
            "Fusion_sequence": "fusion_sequence",
        }
        return cls._process_fusion_caller_rows(path, column_rename, "\t")


class Arriba(FusionCaller):
    """Define parameters for Arriba model"""

    type: Literal[Caller.ARRIBA] = Caller.ARRIBA
    gene1: str = Field(..., description="The 5' gene fusion partner")
    gene2: str = Field(..., description="The 3' gene fusion partner")
    strand1: str = Field(
        ..., description="The strand information for the 5' gene fusion partner"
    )
    strand2: str = Field(
        ..., description="The strand information for the 3' gene fusion partner"
    )
    breakpoint1: str = Field(..., description="The chromosome and breakpoint for gene1")
    breakpoint2: str = Field(..., description="The chromosome and breakpoint for gene2")
    event_type: str = Field(
        ..., description=" An inference about the type of fusion event"
    )
    confidence: str = Field(
        ..., description="A metric describing the confidence of the fusion prediction"
    )
    direction1: str = Field(
        ...,
        description="A description that indicates if the transcript segment starts or ends at breakpoint1",
    )
    direction2: str = Field(
        ...,
        description="A description that indicates if the transcript segment starts or ends at breakpoint2",
    )
    rf: str = Field(
        ...,
        description="A description if the reading frame is preserved for the fusion",
    )
    split_reads1: int = Field(
        ..., description="Number of supporting split fragments with anchor in gene1"
    )
    split_reads2: int = Field(
        ..., description="Number of supporting split fragments with anchor in gene2"
    )
    discordant_mates: int = Field(
        ..., description="Number of discordant mates supporting the fusion"
    )
    coverage1: int = Field(
        ..., description="Number of fragments retained near breakpoint1"
    )
    coverage2: int = Field(
        ..., description="Number of fragments retained near breakpoint2"
    )
    fusion_transcript: str = Field(..., description="The assembled fusion transcript")

    @classmethod
    def load_records(cls, path: Path) -> list["Arriba"]:
        """Load fusions from Arriba tsv file

        :param path: The path to the file of Arriba fusions
        :return A list of Arriba objects, or None if the specified file does not exist
        """
        column_rename = {
            "#gene1": "gene1",
            "strand1(gene/fusion)": "strand1",
            "strand2(gene/fusion)": "strand2",
            "type": "event_type",
            "reading_frame": "rf",
        }
        return cls._process_fusion_caller_rows(path, column_rename, "\t")


class Cicero(FusionCaller):
    """Define parameters for CICERO model"""

    type: Literal[Caller.CICERO] = Caller.CICERO
    gene_5prime: str = Field(..., description="The gene symbol for the 5' partner")
    gene_3prime: str = Field(..., description="The gene symbol for the 3' partner")
    chr_5prime: str = Field(..., description="The chromosome for the 5' partner")
    chr_3prime: str = Field(..., description="The chromosome for the 3' partner")
    pos_5prime: int = Field(
        ..., description="The genomic breakpoint for the 5' partner"
    )
    pos_3prime: int = Field(
        ..., description="The genomic breakpoint for the 3' partner"
    )
    sv_ort: str = Field(
        ...,
        description="Whether the mapping orientation of assembled contig (driven by structural variation) has confident biological meaning",
    )
    event_type: str = Field(
        ...,
        description="The structural variation event that created the called fusion",
    )
    reads_5prime: int = Field(
        ...,
        description="The number of reads that support the breakpoint for the 5' partner",
    )
    reads_3prime: int = Field(
        ...,
        description="The number of reads that support the breakpoint for the 3' partner",
    )
    coverage_5prime: int = Field(
        ..., description="The fragment coverage at the 5' breakpoint"
    )
    coverage_3prime: int = Field(
        ..., description="The fragment coverage at the 3' breakpoint"
    )
    contig: str = Field(..., description="The assembled contig sequence for the fusion")

    @classmethod
    def load_records(cls, path: Path) -> list["Cicero"]:
        """Load fusions from Cicero txt file

        :param path: The path to the file of Cicero fusions
        :return A list of Cicero objects, or None if the specified file does not exist
        """
        column_rename = {
            "geneA": "gene_5prime",
            "geneB": "gene_3prime",
            "chrA": "chr_5prime",
            "chrB": "chr_3prime",
            "posA": "pos_5prime",
            "posB": "pos_3prime",
            "type": "event_type",
            "readsA": "reads_5prime",
            "readsB": "reads_3prime",
            "coverageA": "coverage_5prime",
            "coverageB": "coverage_3prime",
        }
        return cls._process_fusion_caller_rows(path, column_rename, "\t")


class EnFusion(FusionCaller):
    """Define parameters for EnFusion model"""

    type: Literal[Caller.ENFUSION] = Caller.ENFUSION
    gene_5prime: str = Field(..., description="The 5' gene fusion partner")
    gene_3prime: str = Field(..., description="The 3' gene fusion partner")
    chr_5prime: int = Field(..., description="The 5' gene fusion partner chromosome")
    chr_3prime: int = Field(..., description="The 3' gene fusion partner chromosome")
    break_5prime: int = Field(
        ..., description="The 5' gene fusion partner genomic breakpoint"
    )
    break_3prime: int = Field(
        ..., description="The 3' gene fusion partner genomic breakpoint"
    )
    fusion_junction_sequence: str | None = Field(
        None, description="The sequence near the fusion junction"
    )

    @classmethod
    def load_records(cls, path: Path) -> list["EnFusion"]:
        """Load fusions from EnFusion tsv file

        :param path: The path to the file of Enfusion fusions
        :return A list of Enfusion objects, or None if the specified file does not exist
        """
        column_rename = {
            "Gene1": "gene_5prime",
            "Gene2": "gene_3prime",
            "Chr1": "chr_5prime",
            "Chr2": "chr_3prime",
            "Break1": "break_5prime",
            "Break2": "break_3prime",
            "FusionJunctionSequence": "fusion_junction_sequence",
        }
        return cls._process_fusion_caller_rows(path, column_rename, "\t")


class Genie(FusionCaller):
    """Define parameters for Genie model"""

    type: Literal[Caller.GENIE] = Caller.GENIE
    site1_hugo: str = Field(..., description="The HUGO symbol reported at site 1")
    site2_hugo: str = Field(..., description="The HUGO symbol reported at site 2")
    site1_chrom: int = Field(..., description="The chromosome reported at site 1")
    site2_chrom: int = Field(..., description="The chromosome reported at site 2")
    site1_pos: int = Field(..., description="The breakpoint reported at site 1")
    site2_pos: int = Field(..., description="The breakpoint reported at site 2")
    annot: str = Field(..., description="The annotation for the fusion event")
    reading_frame: str = Field(
        ..., description="The reading frame status of the fusion"
    )

    @classmethod
    def load_records(cls, path: Path) -> list["Genie"]:
        """Load fusions from Genie txt file

        :param path: The path to the file of Genie structural variants
        :return A list of Genie objects, or None if the specified file does not exist
        """
        column_rename = {
            "Site1_Hugo_Symbol": "site1_hugo",
            "Site2_Hugo_Symbol": "site2_hugo",
            "Site1_Chromosome": "site1_chrom",
            "Site2_Chromosome": "site2_chrom",
            "Site1_Position": "site1_pos",
            "Site2_Position": "site2_pos",
            "Site2_Effect_On_Frame": "reading_frame",
            "Annotation": "annot",
        }
        return cls._process_fusion_caller_rows(path, column_rename, "\t")
