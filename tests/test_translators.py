"""Module for testing FUSOR Translators"""

import pandas as pd
import pytest

from fusor.models import AssayedFusion


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
        "structural_elements": [
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_152263.4",
                "exon_start": None,
                "exon_start_offset": None,
                "exon_end": 4,
                "exon_end_offset": 5,
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
                            "start": {"type": "Number", "value": 154173078},
                            "end": {"type": "Number", "value": 154173079},
                        },
                    },
                },
            },
            {
                "type": "TranscriptSegmentElement",
                "transcript": "refseq:NM_002609.4",
                "exon_start": 10,
                "exon_start_offset": -559,
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
                            "start": {"type": "Number", "value": 150130527},
                            "end": {"type": "Number", "value": 150130528},
                        },
                    },
                },
                "element_genomic_end": None,
            },
        ],
        "causativeEvent": {"type": "CausativeEvent", "eventType": "rearrangement"},
        "r_frame_preserved": True,
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


def test_gene_element_arriba(translator_instance):
    """Test gene selection for Arriba"""
    genes = "RP1-222H5.1(151985),MIR3672(13973)"
    gene = translator_instance._get_gene_element(genelist=genes, caller="arriba")
    assert gene[0].gene_descriptor.label == "MIR3672"


@pytest.mark.asyncio()
async def test_jaffa(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test JAFFA translator"""
    # Test exonic breakpoint
    jaffa_data = pd.DataFrame(
        [
            [
                "TPM3:PDGFRB",
                "chr1",
                "154170465",
                "-",
                "chr5",
                "150126612",
                "-",
                "TRUE",
                "HighConfidence",
                "TRUE",
            ]
        ],
        columns=[
            "fusion genes",
            "chrom1",
            "base1",
            "strand1",
            "chrom2",
            "base2",
            "strand2",
            "rearrangement",
            "classification",
            "inframe",
        ],
    )
    jaffa_fusor = (await translator_instance.from_jaffa(jaffa_data.iloc[0]))[0].dict()
    # compare_fusions(jaffa_fusor, fusion_data_example.dict())
    assert jaffa_fusor == fusion_data_example.dict()

    # Event Description can be free-string so include seperation assertion statement
    assert jaffa_fusor["causative_event"]["event_description"] == "HighConfidence"
    assert (
        jaffa_fusor["r_frame_preserved"]
        == fusion_data_example.dict()["r_frame_preserved"]
    )

    # Test non-exonic breakpoint
    jaffa_data_nonexonic = pd.DataFrame(
        [
            [
                "TPM3:PDGFRB",
                "chr1",
                "154173079",
                "-",
                "chr5",
                "150130528",
                "-",
                "TRUE",
                "HighConfidence",
                "TRUE",
            ]
        ],
        columns=[
            "fusion genes",
            "chrom1",
            "base1",
            "strand1",
            "chrom2",
            "base2",
            "strand2",
            "rearrangement",
            "classification",
            "inframe",
        ],
    )
    jaffa_fusor_nonexonic = (
        await translator_instance.from_jaffa(jaffa_data_nonexonic.iloc[0])
    )[0].dict()
    compare_fusions(jaffa_fusor_nonexonic, fusion_data_example_nonexonic.dict())

    # Event Description can be free-string so include seperation assertion statement
    assert (
        jaffa_fusor_nonexonic["causative_event"]["event_description"]
        == "HighConfidence"
    )
    assert (
        jaffa_fusor_nonexonic["r_frame_preserved"]
        == fusion_data_example_nonexonic.dict()["r_frame_preserved"]
    )


@pytest.mark.asyncio()
async def test_star_fusion(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test STAR-Fusion translator"""
    # Test exonic breakpoints
    star_fusion_data = pd.DataFrame(
        [
            [
                "TPM3^ENSG00000143549.19",
                "PDGFRB^ENSG00000113721",
                "chr1:154170465:-",
                "chr5:150126612:-",
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

    # Test non-exonic breakpoints
    star_fusion_data_nonexonic = pd.DataFrame(
        [
            [
                "TPM3^ENSG00000143549.19",
                "PDGFRB^ENSG00000113721",
                "chr1:154173079:-",
                "chr5:150130528:-",
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
    star_fusion_fusor_nonexonic = (
        await translator_instance.from_star_fusion(star_fusion_data_nonexonic.iloc[0])
    )[0].dict()
    compare_fusions(star_fusion_fusor_nonexonic, fusion_data_example_nonexonic.dict())
    assert (
        star_fusion_fusor_nonexonic["causative_event"]["event_description"]
        == "INTERCHROMOSOMAL[chr1--chr5]"
    )


@pytest.mark.asyncio()
async def test_fusion_catcher(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test Fusion Catcher translator"""
    # Test exonic breakpoint
    fusion_catcher_data = pd.DataFrame(
        [
            [
                "TPM3",
                "PDGFRB",
                "1:154170465:-",
                "5:150126612:-",
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

    # Test non-exonic breakpoint
    fusion_catcher_data_nonexonic = pd.DataFrame(
        [
            [
                "TPM3",
                "PDGFRB",
                "1:154173079:-",
                "5:150130528:-",
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
    fusion_catcher_fusor_nonexonic = (
        await translator_instance.from_fusion_catcher(
            fusion_catcher_data_nonexonic.iloc[0]
        )
    )[0].dict()
    compare_fusions(
        fusion_catcher_fusor_nonexonic, fusion_data_example_nonexonic.dict()
    )
    assert (
        fusion_catcher_fusor_nonexonic["causative_event"]["event_description"]
        == "exonic(no-known-CDS)/exonic(no-known-CDS)"
    )


@pytest.mark.asyncio()
async def test_fusion_map(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test Fusion Map translator"""
    # Test exonic breakpoint
    fusion_map_data = pd.DataFrame(
        [
            [
                "TPM3",
                "PDGFRB",
                "1",
                "154170465",
                "5",
                "150126612",
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
    assert (
        fusion_map_fusor["r_frame_preserved"]
        == fusion_data_example.dict()["r_frame_preserved"]
    )

    # Test non-exonic breakpoint
    fusion_map_data_nonexonic = pd.DataFrame(
        [
            [
                "TPM3",
                "PDGFRB",
                "1",
                "154173079",
                "5",
                "150130528",
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
    fusion_map_fusor_nonexonic = (
        await translator_instance.from_fusion_map(fusion_map_data_nonexonic.iloc[0])
    )[0].dict()
    compare_fusions(fusion_map_fusor_nonexonic, fusion_data_example_nonexonic.dict())
    assert (
        fusion_map_fusor_nonexonic["causative_event"]["event_description"]
        == "TPM3->PDGFRB,CanonicalPattern[Major],InFrame"
    )
    assert (
        fusion_map_fusor_nonexonic["r_frame_preserved"]
        == fusion_data_example_nonexonic.dict()["r_frame_preserved"]
    )


@pytest.mark.asyncio()
async def test_arriba(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test Arriba translator"""
    # Test exonic breakpoint
    arriba_data = pd.DataFrame(
        [
            [
                "TPM3",
                "PDGFRB",
                "1:154170465",
                "5:150126612",
                "-/-",
                "-/-",
                ".",
                "high",
                "in-frame",
            ]
        ],
        columns=[
            "#gene1",
            "gene2",
            "breakpoint1",
            "breakpoint2",
            "strand1(gene/fusion)",
            "strand2(gene/fusion)",
            "type",
            "confidence",
            "reading_frame",
        ],
    )
    arriba_fusor = (await translator_instance.from_arriba(arriba_data.iloc[0]))[
        0
    ].dict()
    compare_fusions(arriba_fusor, fusion_data_example.dict())
    assert arriba_fusor["causative_event"]["event_description"] == "high"
    assert (
        arriba_fusor["r_frame_preserved"]
        == fusion_data_example.dict()["r_frame_preserved"]
    )

    # Test non-exonic breakpoint
    arriba_data_nonexonic = pd.DataFrame(
        [
            [
                "TPM3",
                "PDGFRB",
                "1:154173079",
                "5:150130528",
                "-/-",
                "-/-",
                ".",
                "high",
                "in-frame",
            ]
        ],
        columns=[
            "#gene1",
            "gene2",
            "breakpoint1",
            "breakpoint2",
            "strand1(gene/fusion)",
            "strand2(gene/fusion)",
            "type",
            "confidence",
            "reading_frame",
        ],
    )
    arriba_fusor_nonexonic = (
        await translator_instance.from_arriba(arriba_data_nonexonic.iloc[0])
    )[0].dict()
    compare_fusions(arriba_fusor_nonexonic, fusion_data_example_nonexonic.dict())
    assert arriba_fusor_nonexonic["causative_event"]["event_description"] == "high"
    assert (
        arriba_fusor_nonexonic["r_frame_preserved"]
        == fusion_data_example_nonexonic.dict()["r_frame_preserved"]
    )


@pytest.mark.asyncio()
async def test_cicero(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test CICERO translator"""
    # Test exonic breakpoint
    cicero_data = pd.DataFrame(
        [["TPM3", "PDGFRB", "1", "-", "154170465", "5", "-", "150126612", "CTX", "1"]],
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
            "frame",
        ],
    )
    cicero_fusor = (await translator_instance.from_cicero(cicero_data.iloc[0]))[
        0
    ].dict()
    compare_fusions(cicero_fusor, fusion_data_example.dict())
    assert cicero_fusor["causative_event"]["event_description"] == "CTX"
    assert (
        cicero_fusor["r_frame_preserved"]
        == fusion_data_example.dict()["r_frame_preserved"]
    )

    # Test non-exonic breakpoint
    cicero_data_nonexonic = pd.DataFrame(
        [["TPM3", "PDGFRB", "1", "-", "154173079", "5", "-", "150130528", "CTX", "1"]],
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
            "frame",
        ],
    )
    cicero_fusor_nonexonic = (
        await translator_instance.from_cicero(cicero_data_nonexonic.iloc[0])
    )[0].dict()
    compare_fusions(cicero_fusor_nonexonic, fusion_data_example_nonexonic.dict())
    assert cicero_fusor["causative_event"]["event_description"] == "CTX"
    assert (
        cicero_fusor_nonexonic["r_frame_preserved"]
        == fusion_data_example_nonexonic.dict()["r_frame_preserved"]
    )


@pytest.mark.asyncio()
async def test_enfusion(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test EnFusion translator"""
    # Test exonic breakpoint
    enfusion_data = pd.DataFrame(
        [["TPM3", "PDGFRB", "1", "-", "154170465", "5", "-", "150126612"]],
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

    # Test non-exonic breakpoint
    enfusion_data_nonexonic = pd.DataFrame(
        [["TPM3", "PDGFRB", "1", "-", "154173079", "5", "-", "150130528"]],
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
    enfusion_fusor_nonexonic = (
        await translator_instance.from_enfusion(enfusion_data_nonexonic.iloc[0])
    )[0].dict()
    compare_fusions(enfusion_fusor_nonexonic, fusion_data_example_nonexonic.dict())


@pytest.mark.asyncio()
async def test_genie(
    fusion_data_example, fusion_data_example_nonexonic, translator_instance
):
    """Test GENIE Translator"""
    # Test exonic breakpoint
    genie_data = pd.DataFrame(
        [
            [
                "TPM3",
                "PDGFRB",
                "1",
                "154170465",
                "5",
                "150126612",
                "exon of TPM3(-): 68 bp before exon 9",
                "exon of PDGFRB(-):1 bp after start of exon 11",
                "TMP3 (NM_152263.4) - PDGFRB (NM_002609.4) fusion",
                "In_frame",
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
            "Site2_Effect_on_Frame",
        ],
    )
    genie_fusor = (await translator_instance.from_genie(genie_data.iloc[0]))[0].dict()
    compare_fusions(genie_fusor, fusion_data_example.dict())
    assert (
        genie_fusor["causative_event"]["event_description"]
        == "TMP3 (NM_152263.4) - PDGFRB (NM_002609.4) fusion"
    )
    assert (
        genie_fusor["r_frame_preserved"]
        == fusion_data_example.dict()["r_frame_preserved"]
    )

    # Test non-exonic breakpoint
    genie_data_nonexonic = pd.DataFrame(
        [
            [
                "TPM3",
                "PDGFRB",
                "1",
                "154173079",
                "5",
                "150130528",
                "exon of TPM3(-): 5 bp after the start of exon 4",
                "exon of PDGFRB(-): 559 bp before the end of exon 10",
                "TMP3 (NM_152263.4) - PDGFRB (NM_002609.4) fusion",
                "In_frame",
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
            "Site2_Effect_on_Frame",
        ],
    )
    genie_fusor_nonexonic = (
        await translator_instance.from_genie(genie_data_nonexonic.iloc[0])
    )[0].dict()
    compare_fusions(genie_fusor_nonexonic, fusion_data_example_nonexonic.dict())
    assert (
        genie_fusor_nonexonic["causative_event"]["event_description"]
        == "TMP3 (NM_152263.4) - PDGFRB (NM_002609.4) fusion"
    )
    assert (
        genie_fusor_nonexonic["r_frame_preserved"]
        == fusion_data_example_nonexonic.dict()["r_frame_preserved"]
    )
