"""Module for translating output from fusion detection methods to FUSOR AssayedFusion
objects
"""

from enum import Enum

import pandas as pd

from fusor.fusor import FUSOR
from fusor.models import (
    Assay,
    AssayedFusion,
    CausativeEvent,
    EventType,
    GeneElement,
    TranscriptSegmentElement,
)


class ReferenceBuild(Enum):
    """Define supported reference builds"""

    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"


class Translator:
    """Class for translating outputs from different fusion detection algorithms
    to FUSOR AssayedFusion objects
    """

    def __init__(self, fusor: FUSOR) -> None:
        """Initialize Translator class
        :param fusor: A FUSOR instance
        """
        self.fusor = fusor

    def _format_fusion(
        self,
        gene_5prime: GeneElement,
        gene_3prime: GeneElement,
        tr_5prime: TranscriptSegmentElement | None = None,
        tr_3prime: TranscriptSegmentElement | None = None,
        ce: CausativeEvent | None = None,
        rf: bool | None = None,
        assay: Assay | None = None,
    ) -> AssayedFusion:
        """Format classes to create AssayedFusion objects
        :param gene_5prime: 5'prime GeneElement
        :param gene_3prime: 3'prime GeneElement
        :param Optional[tr_5prime]: 5'prime TranscriptSegmentElement
        :param Optional[tr_3prime]: 3'prime TranscriptSegmentElement
        :param Optional[ce]: CausativeEvent
        :param Optional[rf]: A boolean indicating if the reading frame is preserved
        :param Optional[assay]: Assay
        :return AssayedFusion object
        """
        params = {}
        if not tr_5prime[0] and not tr_3prime[0]:
            params["structure"] = [gene_5prime[0], gene_3prime[0]]
        elif tr_5prime[0] and not tr_3prime[0]:
            params["structure"] = [tr_5prime[0], gene_3prime[0]]
        elif not tr_5prime[0] and tr_3prime[0]:
            params["structure"] = [gene_5prime[0], tr_3prime[0]]
        else:
            params["structure"] = [tr_5prime[0], tr_3prime[0]]

        if ce:
            params["causativeEvent"] = ce
        if rf:
            params["readingFramePreserved"] = rf
        if assay:
            params["assay"] = assay
        return AssayedFusion(**params), None

    def _get_causative_event(
        self, chrom1: str, chrom2: str, descr: str | None = None
    ) -> CausativeEvent:
        """Infer Causative Event. Currently restricted to rearrangements
        :param chrom1: The chromosome for the 5' partner
        :param chrom2: The chromosome for the 3' partner
        :param annots: An annotation describing the fusion event. This input is supplied to the eventDescription CausativeEvent attribute.
        :return: A CausativeEvent object if construction is successful
        """
        if "rearrangement" in descr:
            return CausativeEvent(
                eventType=EventType("rearrangement"), eventDescription=descr
            )
        if chrom1 != chrom2:
            return CausativeEvent(eventType=EventType("rearrangement"))
        return None

    def _get_gene_element_unnormalized(self, symbol: str) -> GeneElement:
        """Return GeneElement when gene symbol cannot be normalized
        :param symbol: A gene symbol for a fusion partner
        :return: A GeneElement object
        """
        return (
            GeneElement(
                type="GeneElement",
                gene={
                    "id": f"gene:{symbol}",
                    "label": symbol,
                    "type": "Gene",
                },
            ),
            None,
        )

    def _get_gene_element(self, genelist: str, caller: str) -> GeneElement:
        """Return a GeneElement given an individual/list of gene symbols and a
        fusion detection algorithm
        :param genelist: A gene symbol or list of gene symbols
        :param caller: The examined fusion detection algorithm
        :return A GeneElement object
        """
        if "," not in genelist or caller != "arriba":
            ge = self.fusor.gene_element(gene=genelist)
            return ge if ge[0] else self._get_gene_element_unnormalized(genelist)

        genes = genelist.split(",")
        dists = []
        for gene in genes:
            start, end = gene.rfind("("), gene.rfind(")")
            dists.append(int(gene[start + 1 : end]))
        gene = (
            genes[0].split("(")[0] if dists[0] <= dists[1] else genes[1].split("(")[0]
        )
        ge = self.fusor.gene_element(gene=gene)
        return ge if ge[0] else self._get_gene_element_unnormalized(gene)

    def _get_genomic_ac(self, chrom: str, build: ReferenceBuild) -> str:
        """Return a RefSeq genomic accession given a chromosome and a reference build
        :param chrom: A chromosome number
        :param build: The reference build, either GRCh37 or GRCh38
        :return The corresponding refseq genomic accession
        """
        sr = self.fusor.cool_seq_tool.seqrepo_access
        if build == ReferenceBuild.GRCH37:
            alias_list, errors = sr.translate_identifier(f"GRCh37:{chrom}")
        else:
            alias_list, errors = sr.translate_identifier(f"GRCh38:{chrom}")
        if errors:
            return f"Genomic accession for {chrom} could not be retrieved"
        for alias in alias_list:
            if alias.startswith("refseq:"):
                genomic_ac = alias.split(":")[1]
        return genomic_ac

    async def from_jaffa(
        self, jaffa_row: pd.DataFrame, rb: ReferenceBuild
    ) -> AssayedFusion:
        """Parse JAFFA fusion output to create AssayedFusion object
        :param jaffa_row: A row of JAFFA output
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        genes = jaffa_row["fusion genes"].split(":")
        gene_5prime = self._get_gene_element(genes[0], "jaffa")[0].gene.label
        gene_3prime = self._get_gene_element(genes[1], "jaffa")[0].gene.label

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(jaffa_row.chrom1, rb),
            seg_end_genomic=int(jaffa_row.base1),
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(jaffa_row.chrom2, rb),
            seg_start_genomic=int(jaffa_row.base2),
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        if jaffa_row.rearrangement == "TRUE":
            ce = CausativeEvent(
                eventType=EventType("rearrangement"),
                eventDescription=jaffa_row["classification"],
            )
        else:
            ce = None

        rf = bool(jaffa_row.inframe == "TRUE")
        return self._format_fusion(
            gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce, rf
        )

    async def from_star_fusion(
        self, sf_row: pd.DataFrame, rb: ReferenceBuild
    ) -> AssayedFusion:
        """Parse STAR-Fusion output to create AssayedFusion object
        :param sf_row: A row of STAR-Fusion output
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = sf_row["LeftGene"].split("^")[0]
        gene2 = sf_row["RightGene"].split("^")[0]
        gene_5prime = self._get_gene_element(gene1, "star_fusion")[0].gene.label
        gene_3prime = self._get_gene_element(gene2, "star_fusion")[0].gene.label

        five_prime = sf_row["LeftBreakpoint"].split(":")
        three_prime = sf_row["RightBreakpoint"].split(":")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(five_prime[0], rb),
            seg_end_genomic=int(five_prime[1]),
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(three_prime[0], rb),
            seg_start_genomic=int(three_prime[1]),
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        ce = self._get_causative_event(
            five_prime[0], three_prime[0], ",".join(sf_row["annots"])
        )
        return self._format_fusion(gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce)

    async def from_fusion_catcher(
        self, fc_row: pd.DataFrame, rb: ReferenceBuild
    ) -> AssayedFusion:
        """Parse FusionCatcher output to create AssayedFusion object
        :param fc_row: A row of FusionCatcher output
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = fc_row["Gene_1_symbol(5end_fusion_partner)"]
        gene2 = fc_row["Gene_2_symbol(3end_fusion_partner)"]
        gene_5prime = self._get_gene_element(gene1, "fusion_catcher")[0].gene.label
        gene_3prime = self._get_gene_element(gene2, "fusion_catcher")[0].gene.label

        five_prime = fc_row["Fusion_point_for_gene_1(5end_fusion_partner)"].split(":")
        three_prime = fc_row["Fusion_point_for_gene_2(3end_fusion_partner)"].split(":")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(five_prime[0], rb),
            seg_end_genomic=int(five_prime[1]),
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(three_prime[0], rb),
            seg_start_genomic=int(three_prime[1]),
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        ce = self._get_causative_event(
            five_prime[0], three_prime[0], fc_row["Predicted_effect"]
        )
        return self._format_fusion(gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce)

    async def from_fusion_map(
        self, fmap_row: pd.DataFrame, rb: ReferenceBuild
    ) -> AssayedFusion:
        """Parse FusionMap output to create FUSOR AssayedFusion object
        :param fmap_row: A row of FusionMap output
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = fmap_row["KnownGene1"]
        gene2 = fmap_row["KnownGene2"]
        gene_5prime = self._get_gene_element(gene1, "fusion_map")[0].gene.label
        gene_3prime = self._get_gene_element(gene2, "fusion_map")[0].gene.label

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(fmap_row["Chromosome1"], rb),
            seg_end_genomic=int(fmap_row["Position1"]),
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(fmap_row["Chromosome2"], rb),
            seg_start_genomic=int(fmap_row["Position2"]),
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        # Combine columns to create fusion annotation string"
        descr = (
            fmap_row["FusionGene"]
            + ","
            + fmap_row["SplicePatternClass"]
            + ","
            + fmap_row["FrameShiftClass"]
        )
        ce = self._get_causative_event(
            fmap_row["Chromosome1"], fmap_row["Chromosome2"], descr
        )
        rf = bool(fmap_row["FrameShiftClass"] == "InFrame")
        return self._format_fusion(
            gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce, rf
        )

    async def from_arriba(
        self, arriba_row: pd.DataFrame, rb: ReferenceBuild
    ) -> AssayedFusion:
        """Parse Arriba output to create AssayedFusion object
        :param arriba_row: A row of Arriba output
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = arriba_row["#gene1"]
        gene2 = arriba_row["gene2"]

        # Arriba reports two gene symbols if a breakpoint occurs in an intergenic
        # space. We select the gene symbol with the smallest distance from the
        # breakpoint.
        gene_5prime = self._get_gene_element(gene1, "arriba")[0].gene.label
        gene_3prime = self._get_gene_element(gene2, "arriba")[0].gene.label

        breakpoint1 = arriba_row["breakpoint1"].split(":")
        breakpoint2 = arriba_row["breakpoint2"].split(":")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(breakpoint1[0], rb),
            seg_end_genomic=int(breakpoint1[1]),
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(breakpoint2[0], rb),
            seg_start_genomic=int(breakpoint2[1]),
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        ce = (
            CausativeEvent(
                eventType=EventType("read-through"),
                eventDescription=arriba_row["confidence"],
            )
            if "read_through" in arriba_row["type"]
            else CausativeEvent(
                eventType=EventType("rearrangement"),
                eventDescription=arriba_row["confidence"],
            )
        )
        rf = bool(arriba_row["reading_frame"] == "in-frame")
        return self._format_fusion(
            gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce, rf
        )

    async def from_cicero(
        self, cicero_row: pd.DataFrame, rb: ReferenceBuild
    ) -> AssayedFusion:
        """Parse CICERO output to create AssayedFusion object
        :param cicero_row: A row of CICERO output
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = cicero_row["geneA"]
        gene2 = cicero_row["geneB"]
        gene_5prime = self._get_gene_element(gene1, "cicero")[0].gene.label
        gene_3prime = self._get_gene_element(gene2, "cicero")[0].gene.label

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(cicero_row["chrA"], rb),
            seg_end_genomic=int(cicero_row["posA"]),
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(cicero_row["chrB"], rb),
            seg_start_genomic=int(cicero_row["posB"]),
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        if cicero_row["type"] == "read_through":
            ce = CausativeEvent(
                eventType=EventType("read-through"),
                eventDescription=cicero_row["type"],
            )
        else:
            ce = CausativeEvent(
                eventType=EventType("rearrangement"),
                eventDescription=cicero_row["type"],
            )
        rf = bool(int(cicero_row["frame"]) == 1)
        return self._format_fusion(
            gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce, rf
        )

    async def from_mapsplice(
        self, mapsplice_row: pd.DataFrame, rb: ReferenceBuild
    ) -> AssayedFusion:
        """Parse MapSplice output to create AssayedFusion object
        :param mapsplice_row: A row of MapSplice output
        :param rb: The reference build used to call the fusion
        :retun: An AssayedFusion object, if construction is successful
        """
        gene1 = mapsplice_row[60].strip(",")
        gene2 = mapsplice_row[61].strip(",")
        gene_5prime = self._get_gene_element(gene1, "mapsplice")[0].gene.label
        gene_3prime = self._get_gene_element(gene2, "mapsplice")[0].gene.label

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(mapsplice_row[0].split("~")[0], rb),
            seg_end_genomic=int(mapsplice_row[1]),
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(mapsplice_row[0].split("~")[1], rb),
            seg_start_genomic=int(mapsplice_row[2]),
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        ce = self._get_causative_event(
            mapsplice_row[0].split("~")[0], mapsplice_row[0].split("~")[1]
        )
        return self._format_fusion(gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce)

    async def from_enfusion(
        self, enfusion_row: pd.DataFrame, rb: ReferenceBuild
    ) -> AssayedFusion:
        """Parse EnFusion output to create AssayedFusion object
        :param enfusion_row: A row of EnFusion output
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = enfusion_row["Gene1"]
        gene2 = enfusion_row["Gene2"]
        gene_5prime = self._get_gene_element(gene1, "enfusion")[0].gene.label
        gene_3prime = self._get_gene_element(gene2, "enfusion")[0].gene.label

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(enfusion_row["Chr1"], rb),
            seg_end_genomic=int(enfusion_row["Break1"]),
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(enfusion_row["Chr2"], rb),
            seg_start_genomic=int(enfusion_row["Break2"]),
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        ce = self._get_causative_event(enfusion_row["Chr1"], enfusion_row["Chr2"])
        return self._format_fusion(gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce)

    async def from_genie(
        self, genie_row: pd.DataFrame, rb: ReferenceBuild
    ) -> AssayedFusion:
        """Parse GENIE output to create AssayedFusion object
        :param genie_row: A row of EnFusion output
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = genie_row["Site1_Hugo_Symbol"]
        gene2 = genie_row["Site2_Hugo_Symbol"]
        gene_5prime = self._get_gene_element(gene1, "genie")[0].gene.label
        gene_3prime = self._get_gene_element(gene2, "genie")[0].gene.label

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(genie_row["Site1_Chromosome"], rb),
            seg_end_genomic=int(genie_row["Site1_Position"]),
            gene=gene_5prime,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(genie_row["Site2_Chromosome"], rb),
            seg_start_genomic=int(genie_row["Site2_Position"]),
            gene=gene_3prime,
            get_nearest_transcript_junction=True,
        )

        ce = self._get_causative_event(
            genie_row["Site1_Chromosome"],
            genie_row["Site2_Chromosome"],
            genie_row["Annotation"],
        )
        rf = bool(genie_row["Site2_Effect_on_Frame"] == "in frame")
        return self._format_fusion(
            gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce, rf
        )
