"""Module for extracting data from fusion caller output and coverting to pydantic
objects
"""

import csv
import logging
from pathlib import Path

from fusor.fusion_caller_models import (
    JAFFA,
    Arriba,
    Cicero,
    EnFusion,
    FusionCatcher,
    Genie,
    STARFusion,
)

_logger = logging.getLogger(__name__)


def _check_if_file_exists(path: Path) -> bool:
    """Check if fusions file exists

    :param path: The path to the file
    :return ``True`` if the file exists, ``False`` if not
    """
    if not path.exists():
        statement = f"{path!s} does not exist"
        _logger.error(statement)
        return False
    return True


def get_jaffa_records(path: Path) -> list[JAFFA] | None:
    """Load fusions from JAFFA csv file

    :param path: The path to the file of JAFFA fusions
    :return A list of JAFFA objects, or None if the specified file does not exist
    """
    if not _check_if_file_exists(path):
        return None
    fusions_list: list[JAFFA] = []
    column_rename = {
        "fusion genes": "fusion_genes",
        "spanning reads": "spanning_reads",
        "spanning pairs": "spanning_pairs",
    }
    with path.open() as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            row = {column_rename.get(key, key): value for key, value in row.items()}
            fusions_list.append(JAFFA(**row))
    return fusions_list


def get_star_fusion_records(path: Path) -> list[STARFusion] | None:
    """Load fusions from STAR-Fusion tsv file

    :param path: The path to the file of STAR-Fusion fusions
    :return A list of STAR-Fusion objects, or None if the specified file does not exist
    """
    if not _check_if_file_exists(path):
        return None
    fusions_list: list[STARFusion] = []
    column_rename = {
        "LeftGene": "left_gene",
        "RightGene": "right_gene",
        "LeftBreakpoint": "left_breakpoint",
        "RightBreakpoint": "right_breakpoint",
        "JunctionReadCount": "junction_read_count",
        "SpanningFragCount": "spanning_frag_count",
    }
    with path.open() as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            row = {column_rename.get(key, key): value for key, value in row.items()}
            fusions_list.append(STARFusion(**row))
    return fusions_list


def get_fusion_catcher_records(path: Path) -> list[FusionCatcher] | None:
    """Load fusions from FusionCatcher txt file

    :param path: The path to the file of FusionCatcher fusions
    :return A list of FusionCatcher objects, or None if the specified file does not exist
    """
    if not _check_if_file_exists(path):
        return None
    fusions_list: list[FusionCatcher] = []
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
    with path.open() as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            row = {column_rename.get(key, key): value for key, value in row.items()}
            fusions_list.append(FusionCatcher(**row))
    return fusions_list


def get_arriba_records(path: Path) -> list[Arriba] | None:
    """Load fusions from Arriba tsv file

    :param path: The path to the file of Arriba fusions
    :return A list of Arriba objects, or None if the specified file does not exist
    """
    if not _check_if_file_exists(path):
        return None
    fusions_list: list[Arriba] = []
    column_rename = {
        "#gene1": "gene1",
        "strand1(gene/fusion)": "strand1",
        "strand2(gene/fusion)": "strand2",
        "type": "event_type",
        "reading_frame": "rf",
    }
    with path.open() as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            row = {column_rename.get(key, key): value for key, value in row.items()}
            fusions_list.append(Arriba(**row))
    return fusions_list


def get_cicero_records(path: Path) -> list[Cicero] | None:
    """Load fusions from Cicero txt file

    :param path: The path to the file of Cicero fusions
    :return A list of Cicero objects, or None if the specified file does not exist
    """
    if not _check_if_file_exists(path):
        return None
    fusions_list: list[Cicero] = []
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
    with path.open() as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            row = {column_rename.get(key, key): value for key, value in row.items()}
            fusions_list.append(Cicero(**row))
    return fusions_list


def get_enfusion_records(path: Path) -> list[EnFusion] | None:
    """Load fusions from EnFusion tsv file

    :param path: The path to the file of Enfusion fusions
    :return A list of Enfusion objects, or None if the specified file does not exist
    """
    if not _check_if_file_exists(path):
        return None
    fusions_list: list[EnFusion] = []
    column_rename = {
        "Gene1": "gene_5prime",
        "Gene2": "gene_3prime",
        "Chr1": "chr_5prime",
        "Chr2": "chr_3prime",
        "Break1": "break_5prime",
        "Break2": "break_3prime",
        "FusionJunctionSequence": "fusion_junction_sequence",
    }
    with path.open() as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            row = {column_rename.get(key, key): value for key, value in row.items()}
            fusions_list.append(EnFusion(**row))
    return fusions_list


def get_genie_records(path: Path) -> list[Genie] | None:
    """Load fusions from Genie txt file

    :param path: The path to the file of Genie structural variants
    :return A list of Genie objects, or None if the specified file does not exist
    """
    if not _check_if_file_exists(path):
        return None
    fusions_list: list[Genie] = []
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
    with path.open() as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            row = {column_rename.get(key, key): value for key, value in row.items()}
            fusions_list.append(Genie(**row))
    return fusions_list
