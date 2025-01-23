"""Schemas for outputs provided by different fusion callers"""

import csv
from abc import ABC, abstractmethod
from enum import Enum
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, Field


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

    class Config:
        """Allow extra fields from fusion callers to be provided"""

        extra = "allow"

    def _does_file_exist(self, path: Path) -> None:
        """Check if fusions file exists

        :param path: The path to the file
        :return None
        :raise ValueError if the file does not exist at the specified path
        """
        if not path.exists():
            statement = f"{path!s} does not exist"
            raise ValueError(statement)
        return

    def _process_fusion_caller_rows(
        self, path: Path, caller: Caller, column_rename: dict
    ) -> list[Caller]:
        """Convert rows of fusion caller output to Pydantic classes

        :param path: The path to the fusions file
        :param caller: The name of the fusion caller
        :param column_rename: A dictionary of column mappings
        :return: A list of fusions, represented as Pydantic objects
        """
        self._does_file_exist(path)
        fusions_list = []
        with path.open() as csvfile:
            reader = csv.DictReader(
                csvfile, delimiter="," if caller == Caller.JAFFA else "\t"
            )
            for row in reader:
                row = {column_rename.get(key, key): value for key, value in row.items()}
                fusions_list.append(caller(**row))
        return fusions_list

    @abstractmethod
    def load_records(self, path: Path) -> list[Caller]:
        """Abstract method to load records from a fusion caller file."""


class JAFFA(FusionCaller):
    """Define parameters for JAFFA model"""

    type: Literal[Caller.JAFFA] = Caller.JAFFA
    fusion_genes: str = Field(
        None, description="A string containing the two fusion partners"
    )
    chrom1: str = Field(
        None, description="The chromosome indicated in the chrom1 column"
    )
    base1: int = Field(
        None, description="The genomic position indicated in the base1 column"
    )
    chrom2: str = Field(
        None, description="The chromosome indicated in the chrom2 column"
    )
    base2: int = Field(
        None, description="The genomic position indicated in the base2 column"
    )
    rearrangement: bool = Field(
        None, description=" A boolean indicating if a rearrangement occurred"
    )
    classification: str = Field(
        None, description="The classification associated with the called fusion"
    )
    inframe: bool | str = Field(
        None,
        description="A boolean or string indicating if the fusion occurred in-frame",
    )
    spanning_reads: int = Field(
        None,
        description="The number of detected reads that span the junction between the two transcript. Although described as spanning reads, this aligns with our definition of split reads i.e. reads that have sequence belonging to the two fusion partners",
    )
    spanning_pairs: int = Field(
        None,
        description="The number of detected reads that align entirely on either side of the breakpoint",
    )

    def load_records(self, path: Path) -> list["JAFFA"]:
        """Load fusions from JAFFA csv file

        :param path: The path to the file of JAFFA fusions
        :return A list of JAFFA objects, or None if the specified file does not exist
        """
        column_rename = {
            "fusion genes": "fusion_genes",
            "spanning reads": "spanning_reads",
            "spanning pairs": "spanning_pairs",
        }
        return self._process_fusion_caller_rows(path, JAFFA, column_rename)


class STARFusion(BaseModel):
    """Define parameters for STAR-Fusion model"""

    type: Literal[Caller.STAR_FUSION] = Caller.STAR_FUSION
    left_gene: str = Field(..., description="The gene indicated in the LeftGene column")
    right_gene: str = Field(
        None, description="The gene indicated in the RightGene column"
    )
    left_breakpoint: str = Field(
        None, description="The gene indicated in the LeftBreakpoint column"
    )
    right_breakpoint: str = Field(
        None, description="The gene indicated in the RightBreakpoint column"
    )
    annots: str = Field(None, description="The annotations associated with the fusion")
    junction_read_count: int = Field(
        None,
        description="The number of RNA-seq fragments that split the junction between the two transcript segments (from STAR-Fusion documentation)",
    )
    spanning_frag_count: int = Field(
        None,
        description="The number of RNA-seq fragments that encompass the fusion junction such that one read of the pair aligns to a different gene than the other paired-end read of that fragment (from STAR-Fusion documentation)",
    )

    def load_records(self, path: Path) -> list["STARFusion"]:
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
        return self._process_fusion_caller_rows(path, STARFusion, column_rename)


class FusionCatcher(BaseModel):
    """Define parameters for FusionCatcher model"""

    type: Literal[Caller.FUSION_CATCHER] = Caller.FUSION_CATCHER
    five_prime_partner: str = Field(
        None, description="Gene symbol for the 5' fusion partner"
    )
    three_prime_partner: str = Field(
        None, description="Gene symbol for the 3' fusion partner"
    )
    five_prime_fusion_point: str = Field(
        None,
        description="Chromosomal position for the 5' end of the fusion junction. This coordinate is 1-based",
    )
    three_prime_fusion_point: str = Field(
        None,
        description="Chromosomal position for the 3' end of the fusion junction. This coordinate is 1-based",
    )
    predicted_effect: str = Field(
        None,
        description="The predicted effect of the fusion event, created using annotation from the Ensembl database",
    )
    spanning_unique_reads: int = Field(
        None, description="The number of unique reads that map on the fusion junction"
    )
    spanning_reads: int = Field(
        None, description="The number of paired reads that support the fusion"
    )
    fusion_sequence: str = Field(
        None, description="The inferred sequence around the fusion junction"
    )

    def load_records(self, path: Path) -> list["FusionCatcher"]:
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
        return self._process_fusion_caller_rows(path, FusionCatcher, column_rename)


class Arriba(BaseModel):
    """Define parameters for Arriba model"""

    type: Literal[Caller.ARRIBA] = Caller.ARRIBA
    gene1: str = Field(..., description="The 5' gene fusion partner")
    gene2: str = Field(..., description="The 3' gene fusion partner")
    strand1: str = Field(
        None, description="The strand information for the 5' gene fusion partner"
    )
    strand2: str = Field(
        None, description="The strand information for the 3' gene fusion partner"
    )
    breakpoint1: str = Field(
        None, description="The chromosome and breakpoint for gene1"
    )
    breakpoint2: str = Field(
        None, description="The chromosome and breakpoint for gene2"
    )
    event_type: str = Field(
        None, description=" An inference about the type of fusion event"
    )
    confidence: str = Field(
        None, description="A metric describing the confidence of the fusion prediction"
    )
    direction1: str = Field(
        None,
        description="A description that indicates if the transcript segment starts or ends at breakpoint1",
    )
    direction2: str = Field(
        None,
        description="A description that indicates if the transcript segment starts or ends at breakpoint2",
    )
    rf: str = Field(
        None,
        description="A description if the reading frame is preserved for the fusion",
    )
    split_reads1: int = Field(
        None, description="Number of supporting split fragments with anchor in gene1"
    )
    split_reads2: int = Field(
        None, description="Number of supporting split fragments with anchor in gene2"
    )
    discordant_mates: int = Field(
        None, description="Number of discordant mates supporting the fusion"
    )
    coverage1: int = Field(
        None, description="Number of fragments retained near breakpoint1"
    )
    coverage2: int = Field(
        None, description="Number of fragments retained near breakpoint2"
    )
    fusion_transcript: str = Field(None, description="The assembled fusion transcript")

    def load_records(self, path: Path) -> list["Arriba"]:
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
        return self._process_fusion_caller_rows(path, Arriba, column_rename)


class Cicero(BaseModel):
    """Define parameters for CICERO model"""

    type: Literal[Caller.CICERO] = Caller.CICERO
    gene_5prime: str = Field(None, description="The gene symbol for the 5' partner")
    gene_3prime: str = Field(None, description="The gene symbol for the 3' partner")
    chr_5prime: str = Field(None, description="The chromosome for the 5' partner")
    chr_3prime: str = Field(None, description="The chromosome for the 3' partner")
    pos_5prime: int = Field(
        None, description="The genomic breakpoint for the 5' partner"
    )
    pos_3prime: int = Field(
        None, description="The genomic breakpoint for the 3' partner"
    )
    sv_ort: str = Field(
        None,
        description="Whether the mapping orientation of assembled contig (driven by structural variation) has confident biological meaning",
    )
    event_type: str = Field(
        None,
        description="The structural variation event that created the called fusion",
    )
    reads_5prime: int = Field(
        None,
        description="The number of reads that support the breakpoint for the 5' partner",
    )
    reads_3prime: int = Field(
        None,
        description="The number of reads that support the breakpoint for the 3' partner",
    )
    coverage_5prime: int = Field(
        None, description="The fragment coverage at the 5' breakpoint"
    )
    coverage_3prime: int = Field(
        None, description="The fragment coverage at the 3' breakpoint"
    )
    contig: str = Field(..., description="The assembled contig sequence for the fusion")

    def load_records(self, path: Path) -> list["Cicero"]:
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
        return self._process_fusion_caller_rows(path, Cicero, column_rename)


class EnFusion(BaseModel):
    """Define parameters for EnFusion model"""

    type: Literal[Caller.ENFUSION] = Caller.ENFUSION
    gene_5prime: str = Field(None, description="The 5' gene fusion partner")
    gene_3prime: str = Field(None, description="The 3' gene fusion partner")
    chr_5prime: int = Field(None, description="The 5' gene fusion partner chromosome")
    chr_3prime: int = Field(None, description="The 3' gene fusion partner chromosome")
    break_5prime: int = Field(
        None, description="The 5' gene fusion partner genomic breakpoint"
    )
    break_3prime: int = Field(
        None, description="The 3' gene fusion partner genomic breakpoint"
    )
    fusion_junction_sequence: str | None = Field(
        None, description="The sequence near the fusion junction"
    )

    def load_records(self, path: Path) -> list["EnFusion"]:
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
        return self._process_fusion_caller_rows(path, EnFusion, column_rename)


class Genie(BaseModel):
    """Define parameters for Genie model"""

    type: Literal[Caller.GENIE] = Caller.GENIE
    site1_hugo: str = Field(None, description="The HUGO symbol reported at site 1")
    site2_hugo: str = Field(None, description="The HUGO symbol reported at site 2")
    site1_chrom: int = Field(None, description="The chromosome reported at site 1")
    site2_chrom: int = Field(None, description="The chromosome reported at site 2")
    site1_pos: int = Field(None, description="The breakpoint reported at site 1")
    site2_pos: int = Field(None, description="The breakpoint reported at site 2")
    annot: str = Field(None, description="The annotation for the fusion event")
    reading_frame: str = Field(
        None, description="The reading frame status of the fusion"
    )

    def load_records(self, path: Path) -> list["Genie"]:
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
        return self._process_fusion_caller_rows(path, Genie, column_rename)
