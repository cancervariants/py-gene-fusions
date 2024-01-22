"""Module for testing FUSOR Translators"""
import pandas as pd
import pytest

from fusor.models import AssayedFusion
from fusor.translator import Translator


@pytest.fixture(scope="module")
def fusion_data_example():
    """Create example assayed fusion TPM3::PDGFRB"""
    params = {
        "type": "AssayedFusion",
        "structural_elements": [
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_152263.4",
                "exon_start": None,
                "exon_start_offset": None,
                "exon_end": 8,
                "exon_end_offset": -68,
                "gene_descriptor": {
                    "id": "normalize.gene:TPM3",
                    "type": "GeneDescriptor",
                    "label": "TPM3",
                    "gene_id": "hgnc:12012",
                },
                "element_genomic_start": None,
                "element_genomic_end": {
                    "id": "fusor.location_descriptor:NC_000001.11",
                    "type": "LocationDescriptor",
                    "label": "NC_000001.11",
                    "location": {
                        "type": "SequenceLocation",
                        "sequence_id": "refseq:NC_000001.11",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 154170467},
                            "end": {"type": "Number", "value": 154170468},
                        },
                    },
                },
            },
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_002609.4",
                "exon_start": 11,
                "exon_start_offset": 1,
                "exon_end": None,
                "exon_end_offset": None,
                "gene_descriptor": {
                    "id": "normalize.gene:PDGFRB",
                    "type": "GeneDescriptor",
                    "label": "PDGFRB",
                    "gene_id": "hgnc:8804",
                },
                "element_genomic_start": {
                    "id": "fusor.location_descriptor:NC_000005.10",
                    "type": "LocationDescriptor",
                    "label": "NC_000005.10",
                    "location": {
                        "type": "SequenceLocation",
                        "sequence_id": "refseq:NC_000005.10",
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": 150126613},
                            "end": {"type": "Number", "value": 150126614},
                        },
                    },
                },
                "element_genomic_end": None,
            },
        ],
        "causative_event": {"type": "CausativeEvent", "event_type": "rearrangement"},
        "assay": None,
    }
    return AssayedFusion(**params)


def compare_fusions(alg_fusion: dict, fusion_data_example: dict):
    """Create helper function for comparing fusion output"""
    assert (
        alg_fusion["structural_elements"][0]["transcript"]
        == fusion_data_example["structural_elements"][0]["transcript"]
    )
    assert (
        alg_fusion["structural_elements"][0]["exon_start"]
        == fusion_data_example["structural_elements"][0]["exon_start"]
    )
    assert (
        alg_fusion["structural_elements"][0]["exon_start_offset"]
        == fusion_data_example["structural_elements"][0]["exon_start_offset"]
    )
    assert (
        alg_fusion["structural_elements"][0]["exon_end"]
        == fusion_data_example["structural_elements"][0]["exon_end"]
    )
    assert (
        alg_fusion["structural_elements"][0]["exon_end_offset"]
        == fusion_data_example["structural_elements"][0]["exon_end_offset"]
    )
    assert (
        alg_fusion["structural_elements"][0]["gene_descriptor"]
        == fusion_data_example["structural_elements"][0]["gene_descriptor"]
    )
    assert (
        alg_fusion["structural_elements"][0]["element_genomic_start"]
        == fusion_data_example["structural_elements"][0]["element_genomic_start"]
    )
    assert (
        alg_fusion["structural_elements"][0]["element_genomic_end"]
        == fusion_data_example["structural_elements"][0]["element_genomic_end"]
    )
    assert (
        alg_fusion["structural_elements"][1]["transcript"]
        == fusion_data_example["structural_elements"][1]["transcript"]
    )
    assert (
        alg_fusion["structural_elements"][1]["exon_start"]
        == fusion_data_example["structural_elements"][1]["exon_start"]
    )
    assert (
        alg_fusion["structural_elements"][1]["exon_start_offset"]
        == fusion_data_example["structural_elements"][1]["exon_start_offset"]
    )
    assert (
        alg_fusion["structural_elements"][1]["exon_end"]
        == fusion_data_example["structural_elements"][1]["exon_end"]
    )
    assert (
        alg_fusion["structural_elements"][1]["exon_end_offset"]
        == fusion_data_example["structural_elements"][1]["exon_end_offset"]
    )
    assert (
        alg_fusion["structural_elements"][1]["gene_descriptor"]
        == fusion_data_example["structural_elements"][1]["gene_descriptor"]
    )
    assert (
        alg_fusion["structural_elements"][1]["element_genomic_start"]
        == fusion_data_example["structural_elements"][1]["element_genomic_start"]
    )
    assert (
        alg_fusion["structural_elements"][1]["element_genomic_end"]
        == fusion_data_example["structural_elements"][1]["element_genomic_end"]
    )
    assert (
        alg_fusion["causative_event"]["event_type"]
        == fusion_data_example["causative_event"]["event_type"]
    )
    assert alg_fusion["assay"] == fusion_data_example["assay"]


@pytest.mark.asyncio
async def test_jaffa(fusion_data_example, fusor_instance):
    """Test JAFFA translator"""
    translator_instance = Translator(fusor_instance)
    jaffa_data = pd.DataFrame(
        [["TPM3:PDGFRB", "chr1", "154170400", "-", "chr5", "150125578", "-", True]],
        columns=[
            "fusion genes",
            "chrom1",
            "base1",
            "strand1",
            "chrom2",
            "base2",
            "strand2",
            "rearrangement",
        ],
    )
    jaffa_fusor = (await translator_instance.from_jaffa(jaffa_data.iloc[0]))[0].dict()
    compare_fusions(jaffa_fusor, fusion_data_example.dict())


@pytest.mark.asyncio
async def test_star_fusion(fusion_data_example, fusor_instance):
    """Test STAR-Fusion translator"""
    translator_instance = Translator(fusor_instance)
    star_fusion_data = pd.DataFrame(
        [
            [
                "TPM3^ENSG00000143549.19",
                "PDGFRB^ENSG00000113721",
                "chr1:154170400:-",
                "chr5:150125578:-",
                "INTERCHROMOSOMAL[chr1--chr5]",
            ]
        ],
        columns=[
            "LeftGene",
            "RightGene",
            "LeftBreakpoint",
            "RightBreakpoint",
            "annots",
        ],
    )
    star_fusion_fusor = (
        await translator_instance.from_star_fusion(star_fusion_data.iloc[0])
    )[0].dict()
    compare_fusions(star_fusion_fusor, fusion_data_example.dict())
    assert (
        star_fusion_fusor["causative_event"]["event_description"]
        == "INTERCHROMOSOMAL[chr1--chr5]"
    )


@pytest.mark.asyncio
async def test_fusion_catcher(fusion_data_example, fusor_instance):
    """Test Fusion Catcher translator"""
    translator_instance = Translator(fusor_instance)
    fusion_catcher_data = pd.DataFrame(
        [
            [
                "TPM3",
                "PDGFRB",
                "1:154170400:-",
                "5:150125578:-",
                "exonic(no-known-CDS)/exonic(no-known-CDS)",
            ]
        ],
        columns=[
            "Gene_1_symbol(5end_fusion_partner)",
            "Gene_2_symbol(3end_fusion_partner)",
            "Fusion_point_for_gene_1(5end_fusion_partner)",
            "Fusion_point_for_gene_2(3end_fusion_partner)",
            "Predicted_effect",
        ],
    )
    fusion_catcher_fusor = (
        await translator_instance.from_fusion_catcher(fusion_catcher_data.iloc[0])
    )[0].dict()
    compare_fusions(fusion_catcher_fusor, fusion_data_example.dict())
    assert (
        fusion_catcher_fusor["causative_event"]["event_description"]
        == "exonic(no-known-CDS)/exonic(no-known-CDS)"
    )


@pytest.mark.asyncio
async def test_fusion_map(fusion_data_example, fusor_instance):
    """Test Fusion Map translator"""
    translator_instance = Translator(fusor_instance)
    fusion_map_data = pd.DataFrame(
        [
            [
                "TPM3",
                "PDGFRB",
                "1",
                "154170400",
                "5",
                "150125578",
                "--",
                "TPM3->PDGFRB",
                "CanonicalPattern[Major]",
                "InFrame",
            ]
        ],
        columns=[
            "KnownGene1",
            "KnownGene2",
            "Chromosome1",
            "Position1",
            "Chromosome2",
            "Position2",
            "Strand",
            "FusionGene",
            "SplicePatternClass",
            "FrameShiftClass",
        ],
    )
    fusion_map_fusor = (
        await translator_instance.from_fusion_map(fusion_map_data.iloc[0])
    )[0].dict()
    compare_fusions(fusion_map_fusor, fusion_data_example.dict())
    assert (
        fusion_map_fusor["causative_event"]["event_description"]
        == "TPM3->PDGFRB,CanonicalPattern[Major],InFrame"
    )


@pytest.mark.asyncio
async def test_arriba(fusion_data_example, fusor_instance):
    """Test Arriba translator"""
    translator_instance = Translator(fusor_instance)
    arriba_data = pd.DataFrame(
        [["TPM3", "PDGFRB", "1:154170400", "5:150125578", "-/-", "-/-", "."]],
        columns=[
            "gene1",
            "gene2",
            "breakpoint1",
            "breakpoint2",
            "strand1(gene/fusion)",
            "strand2(gene/fusion)",
            "type",
        ],
    )
    arriba_fusor = (await translator_instance.from_arriba(arriba_data.iloc[0]))[
        0
    ].dict()
    compare_fusions(arriba_fusor, fusion_data_example.dict())


@pytest.mark.asyncio
async def test_cicero(fusion_data_example, fusor_instance):
    """Test CICERO translator"""
    translator_instance = Translator(fusor_instance)
    cicero_data = pd.DataFrame(
        [["TPM3", "PDGFRB", "1", "-", "154170400", "5", "-", "150125578", "CTX"]],
        columns=[
            "geneA",
            "geneB",
            "chrA",
            "ortA",
            "posA",
            "chrB",
            "ortB",
            "posB",
            "type",
        ],
    )
    cicero_fusor = (await translator_instance.from_cicero(cicero_data.iloc[0]))[
        0
    ].dict()
    compare_fusions(cicero_fusor, fusion_data_example.dict())


@pytest.mark.asyncio
async def test_enfusion(fusion_data_example, fusor_instance):
    """Test EnFusion translator"""
    translator_instance = Translator(fusor_instance)
    enfusion_data = pd.DataFrame(
        [["TPM3", "PDGFRB", "1", "-", "154170400", "5", "-", "150125578"]],
        columns=[
            "Gene1",
            "Gene2",
            "Chr1",
            "Strand1",
            "Break1",
            "Chr2",
            "Strand2",
            "Break2",
        ],
    )
    enfusion_fusor = (await translator_instance.from_enfusion(enfusion_data.iloc[0]))[
        0
    ].dict()
    compare_fusions(enfusion_fusor, fusion_data_example.dict())


@pytest.mark.asyncio
async def test_genie(fusion_data_example, fusor_instance):
    """Test GENIE Translator"""
    translator_instance = Translator(fusor_instance)
    genie_data = pd.DataFrame(
        [
            [
                "TPM3",
                "PDGFRB",
                "1",
                "154170400",
                "5",
                "150125578",
                "exon of TPM3(-): 68bp before exon 9",
                "exon of PDGFRB(-):1 bp after start of exon 11",
                "TMP3 (NM_152263.4) - PDGFRB (NM_002609.4) fusion",
            ]
        ],
        columns=[
            "Site1_Hugo_Symbol",
            "Site2_Hugo_Symbol",
            "Site1_Chromosome",
            "Site1_Position",
            "Site2_Chromosome",
            "Site2_Position",
            "Site1_Description",
            "Site2_Description",
            "Annotation",
        ],
    )
    genie_fusor = (await translator_instance.from_genie(genie_data.iloc[0]))[0].dict()
    compare_fusions(genie_fusor, fusion_data_example.dict())
    assert (
        genie_fusor["causative_event"]["event_description"]
        == "TMP3 (NM_152263.4) - PDGFRB (NM_002609.4) fusion"
    )
