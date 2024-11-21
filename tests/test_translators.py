"""Module for testing FUSOR Translators"""

import polars as pl
import pytest
from cool_seq_tool.schemas import Assembly

from fusor.models import AssayedFusion
from fusor.translator import Caller


@pytest.fixture(scope="module")
def fusion_data_example():
    """Create example assayed fusion for TPM3::PDGFRB with exonic breakpoints"""
    params = {
        "type": "AssayedFusion",
        "structure": [
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_152263.4",
                "exonEnd": 8,
                "exonEndOffset": -66,
                "gene": {"id": "hgnc:12012", "type": "Gene", "label": "TPM3"},
                "elementGenomicEnd": {
                    "id": "ga4gh:SL.6lXn5i3zqcZUfmtBSieTiVL4Nt2gPGKY",
                    "type": "SequenceLocation",
                    "digest": "6lXn5i3zqcZUfmtBSieTiVL4Nt2gPGKY",
                    "sequenceReference": {
                        "id": "refseq:NC_000001.11",
                        "type": "SequenceReference",
                        "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                    },
                    "start": 154170465,
                },
            },
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_002609.4",
                "exonStart": 11,
                "exonStartOffset": 2,
                "gene": {"id": "hgnc:8804", "type": "Gene", "label": "PDGFRB"},
                "elementGenomicStart": {
                    "id": "ga4gh:SL.Sp1lwuHbRCkWIoe4zzwVKPsS8zK8i0ck",
                    "type": "SequenceLocation",
                    "digest": "Sp1lwuHbRCkWIoe4zzwVKPsS8zK8i0ck",
                    "sequenceReference": {
                        "id": "refseq:NC_000005.10",
                        "type": "SequenceReference",
                        "refgetAccession": "SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI",
                    },
                    "end": 150126612,
                },
            },
        ],
        "causativeEvent": {"type": "CausativeEvent", "eventType": "rearrangement"},
        "r_frame_preserved": True,
        "assay": None,
    }
    return AssayedFusion(**params)


@pytest.fixture(scope="module")
def fusion_data_example_nonexonic():
    """Create example assayed fusion for TPM3::PDGFRB with non-exonic breakpoints"""
    params = {
        "type": "AssayedFusion",
        "structure": [
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_152263.4",
                "exonEnd": 4,
                "exonEndOffset": 5,
                "gene": {"id": "hgnc:12012", "type": "Gene", "label": "TPM3"},
                "elementGenomicEnd": {
                    "id": "ga4gh:SL.O1rVKQA2FTdy_FFWg3qJVSTG_TF_Mkex",
                    "type": "SequenceLocation",
                    "digest": "O1rVKQA2FTdy_FFWg3qJVSTG_TF_Mkex",
                    "sequenceReference": {
                        "id": "refseq:NC_000001.11",
                        "type": "SequenceReference",
                        "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                    },
                    "start": 154173078,
                },
            },
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_002609.4",
                "exonStart": 11,
                "exonStartOffset": -559,
                "gene": {"id": "hgnc:8804", "type": "Gene", "label": "PDGFRB"},
                "elementGenomicStart": {
                    "id": "ga4gh:SL.GtoWMuox4tOyX2I5L9Baobnpgc1pDIVJ",
                    "type": "SequenceLocation",
                    "digest": "GtoWMuox4tOyX2I5L9Baobnpgc1pDIVJ",
                    "sequenceReference": {
                        "id": "refseq:NC_000005.10",
                        "type": "SequenceReference",
                        "refgetAccession": "SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI",
                    },
                    "end": 150127173,
                },
            },
        ],
        "causativeEvent": {"type": "CausativeEvent", "eventType": "rearrangement"},
        "r_frame_preserved": True,
        "assay": None,
    }
    return AssayedFusion(**params)


def test_gene_element_arriba(translator_instance):
    """Test gene selection for Arriba"""
    genes = "RP1-222H5.1(151985),MIR3672(13973)"
    gene = translator_instance._get_gene_element(genes=genes, caller=Caller.ARRIBA)
    assert gene[0].gene.label == "MIR3672"


@pytest.mark.asyncio()
async def test_jaffa(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test JAFFA translator"""
    # Test exonic breakpoint
    fusion_genes = "TPM3:PDGFRB"
    chrom1 = "chr1"
    base1 = 154170465
    chrom2 = "chr5"
    base2 = 150126612
    rearrangement = True
    classification = "HighConfidence"
    inframe = True

    jaffa_fusor = await translator_instance.from_jaffa(
        fusion_genes,
        chrom1,
        base1,
        chrom2,
        base2,
        rearrangement,
        classification,
        inframe,
        Assembly.GRCH38.value,
    )
    assert jaffa_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    fusion_genes = "TPM3:PDGFRB"
    chrom1 = "chr1"
    base1 = 154173078
    chrom2 = "chr5"
    base2 = 150127173
    rearrangement = True
    classification = "HighConfidence"
    inframe = True

    jaffa_fusor_nonexonic = await translator_instance.from_jaffa(
        fusion_genes,
        chrom1,
        base1,
        chrom2,
        base2,
        rearrangement,
        classification,
        inframe,
        Assembly.GRCH38.value,
    )
    assert jaffa_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure


@pytest.mark.asyncio()
async def test_star_fusion(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test STAR-Fusion translator"""
    # Test exonic breakpoints
    left_gene = "TPM3^ENSG00000143549.19"
    right_gene = "PDGFRB^ENSG00000113721"
    left_breakpoint = "chr1:154170465:-"
    right_breakpoint = "chr5:150126612:-"
    annots = '["INTRACHROMOSOMAL[chr16:0.23Mb]"]'

    star_fusion_fusor = await translator_instance.from_star_fusion(
        left_gene,
        right_gene,
        left_breakpoint,
        right_breakpoint,
        annots,
        Assembly.GRCH38.value,
    )
    assert star_fusion_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoints
    left_gene = "TPM3^ENSG00000143549.19"
    right_gene = "PDGFRB^ENSG00000113721"
    left_breakpoint = "chr1:154173078:-"
    right_breakpoint = "chr5:150127173:-"
    annots = '["INTRACHROMOSOMAL[chr16:0.23Mb]"]'

    star_fusion_fusor_nonexonic = await translator_instance.from_star_fusion(
        left_gene,
        right_gene,
        left_breakpoint,
        right_breakpoint,
        annots,
        Assembly.GRCH38.value,
    )
    assert (
        star_fusion_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure
    )


@pytest.mark.asyncio()
async def test_fusion_catcher(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test Fusion Catcher translator"""
    # Test exonic breakpoint
    five_prime_partner = "TPM3"
    three_prime_partner = "PDGFRB"
    five_prime_fusion_point = "1:154170466:-"
    three_prime_fusion_point = "5:150126613:-"
    predicted_effect = "exonic(no-known-CDS)/exonic(no-known-CDS)"

    fusion_catcher_fusor = await translator_instance.from_fusion_catcher(
        five_prime_partner,
        three_prime_partner,
        five_prime_fusion_point,
        three_prime_fusion_point,
        predicted_effect,
        Assembly.GRCH38.value,
    )
    assert fusion_catcher_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    five_prime_partner = "TPM3"
    three_prime_partner = "PDGFRB"
    five_prime_fusion_point = "1:154173079:-"
    three_prime_fusion_point = "5:150127174:-"
    predicted_effect = "exonic(no-known-CDS)/exonic(no-known-CDS)"

    fusion_catcher_fusor_nonexonic = await translator_instance.from_fusion_catcher(
        five_prime_partner,
        three_prime_partner,
        five_prime_fusion_point,
        three_prime_fusion_point,
        predicted_effect,
        Assembly.GRCH38.value,
    )
    assert (
        fusion_catcher_fusor_nonexonic.structure
        == fusion_data_example_nonexonic.structure
    )


@pytest.mark.asyncio()
async def test_fusion_map(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test Fusion Map translator"""
    # Test exonic breakpoint
    fusion_map_data = pl.DataFrame(
        {
            "KnownGene1": "TPM3",
            "KnownGene2": "PDGFRB",
            "Chromosome1": "1",
            "Position1": "154170465",
            "Chromosome2": "5",
            "Position2": "150126612",
            "FusionGene": "TPM3->PDGFRB",
            "SplicePatternClass": "CanonicalPattern[Major]",
            "FrameShiftClass": "InFrame",
        }
    )
    fusion_map_fusor = await translator_instance.from_fusion_map(
        fusion_map_data, Assembly.GRCH38.value
    )
    assert fusion_map_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    fusion_map_data_nonexonic = pl.DataFrame(
        {
            "KnownGene1": "TPM3",
            "KnownGene2": "PDGFRB",
            "Chromosome1": "1",
            "Position1": "154173078",
            "Chromosome2": "5",
            "Position2": "150127173",
            "FusionGene": "TPM3->PDGFRB",
            "SplicePatternClass": "CanonicalPattern[Major]",
            "FrameShiftClass": "InFrame",
        }
    )
    fusion_map_fusor_nonexonic = await translator_instance.from_fusion_map(
        fusion_map_data_nonexonic, Assembly.GRCH38.value
    )
    assert (
        fusion_map_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure
    )


@pytest.mark.asyncio()
async def test_arriba(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test Arriba translator"""
    # Test exonic breakpoint
    arriba_data = pl.DataFrame(
        {
            "#gene1": "TPM3",
            "gene2": "PDGFRB",
            "breakpoint1": "1:154170465",
            "breakpoint2": "5:150126612",
            "type": ".",
            "confidence": "high",
            "reading_frame": "in-frame",
        }
    )
    arriba_fusor = await translator_instance.from_arriba(
        arriba_data, Assembly.GRCH38.value
    )
    assert arriba_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    arriba_data_nonexonic = pl.DataFrame(
        {
            "#gene1": "TPM3",
            "gene2": "PDGFRB",
            "breakpoint1": "1:154173078",
            "breakpoint2": "5:150127173",
            "type": ".",
            "confidence": "high",
            "reading_frame": "in-frame",
        }
    )
    arriba_fusor_nonexonic = await translator_instance.from_arriba(
        arriba_data_nonexonic, Assembly.GRCH38.value
    )
    assert arriba_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure


@pytest.mark.asyncio()
async def test_cicero(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test CICERO translator"""
    # Test exonic breakpoint
    cicero_data = pl.DataFrame(
        {
            "geneA": "TPM3",
            "geneB": "PDGFRB",
            "chrA": "1",
            "posA": "154170465",
            "chrB": "5",
            "posB": "150126612",
            "type": "CTX",
        }
    )
    cicero_fusor = await translator_instance.from_cicero(
        cicero_data, Assembly.GRCH38.value
    )
    assert cicero_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    cicero_data_nonexonic = pl.DataFrame(
        {
            "geneA": "TPM3",
            "geneB": "PDGFRB",
            "chrA": "1",
            "posA": "154173078",
            "chrB": "5",
            "posB": "150127173",
            "type": "CTX",
        }
    )
    cicero_fusor_nonexonic = await translator_instance.from_cicero(
        cicero_data_nonexonic, Assembly.GRCH38.value
    )
    assert cicero_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure


@pytest.mark.asyncio()
async def test_enfusion(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test EnFusion translator"""
    # Test exonic breakpoint
    enfusion_data = pl.DataFrame(
        {
            "Gene1": "TPM3",
            "Gene2": "PDGFRB",
            "Chr1": "1",
            "Break1": "154170465",
            "Chr2": "5",
            "Break2": "150126612",
        }
    )
    enfusion_fusor = await translator_instance.from_enfusion(
        enfusion_data, Assembly.GRCH38.value
    )
    assert enfusion_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    enfusion_data_nonexonic = pl.DataFrame(
        {
            "Gene1": "TPM3",
            "Gene2": "PDGFRB",
            "Chr1": "1",
            "Break1": "154173078",
            "Chr2": "5",
            "Break2": "150127173",
        }
    )
    enfusion_fusor_nonexonic = await translator_instance.from_enfusion(
        enfusion_data_nonexonic, Assembly.GRCH38.value
    )
    assert enfusion_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure


@pytest.mark.asyncio()
async def test_genie(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test GENIE Translator"""
    # Test exonic breakpoint
    genie_data = pl.DataFrame(
        {
            "Site1_Hugo_Symbol": "TPM3",
            "Site2_Hugo_Symbol": "PDGFRB",
            "Site1_Chromosome": "1",
            "Site1_Position": "154170465",
            "Site2_Chromosome": "5",
            "Site2_Position": "150126612",
            "Annotation": "TMP3 (NM_152263.4) - PDGFRB (NM_002609.4) fusion",
            "Site2_Effect_on_Frame": "In_frame",
        }
    )
    genie_fusor = await translator_instance.from_genie(
        genie_data, Assembly.GRCH38.value
    )
    assert genie_fusor.structure == fusion_data_example.structure

    # Test non-exonic breakpoint
    genie_data_nonexonic = pl.DataFrame(
        {
            "Site1_Hugo_Symbol": "TPM3",
            "Site2_Hugo_Symbol": "PDGFRB",
            "Site1_Chromosome": "1",
            "Site1_Position": "154173078",
            "Site2_Chromosome": "5",
            "Site2_Position": "150127173",
            "Annotation": "TMP3 (NM_152263.4) - PDGFRB (NM_002609.4) fusion",
            "Site2_Effect_on_Frame": "In_frame",
        }
    )
    genie_fusor_nonexonic = await translator_instance.from_genie(
        genie_data_nonexonic, Assembly.GRCH38.value
    )
    assert genie_fusor_nonexonic.structure == fusion_data_example_nonexonic.structure
