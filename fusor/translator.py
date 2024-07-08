"""Module for translating output from fusion detection methods to FUSOR AssayedFusion
objects
"""

from typing import Optional

import pandas as pd
from cool_seq_tool.schemas import Strand

from fusor.fusor import FUSOR
from fusor.models import (
    AssayedFusion,
    CausativeEvent,
    EventType,
    GeneDescriptor,
    GeneElement,
    TranscriptSegmentElement,
)


class Translator:
    """Class for translating outputs from different fusion detection algorithms
    to FUSOR AssayedFusion objects
    """

    def __init__(self, fusor: FUSOR) -> None:
        """Initialize Translator class
        :param fusor: A FUSOR instance
        """
        self.fusor = FUSOR()

    def _format_fusion(
        self,
        gene_5prime: GeneElement,
        gene_3prime: GeneElement,
        tr_5prime: Optional[TranscriptSegmentElement] = None,
        tr_3prime: Optional[TranscriptSegmentElement] = None,
        ce: Optional[CausativeEvent] = None,
        rf: Optional[bool] = None,
    ) -> AssayedFusion:
        """Format classes to create AssayedFusion objects
        :param gene_5prime: 5'prime GeneElement
        :param gene_3prime: 3'prime GeneElement
        :param Optional[tr_5prime]: 5'prime TranscriptSegmentElement
        :param Optional[tr_3prime]: 3'prime TranscriptSegmentElement
        :param Optional[ce]: CausativeEvent
        :param Optional[rf]: A boolean indicating if the reading frame is preserved
        :return AssayedFusion object
        """
        params = {}
        if not tr_5prime[0] and not tr_3prime[0]:
            params["structural_elements"] = [gene_5prime[0], gene_3prime[0]]
        elif tr_5prime[0] and not tr_3prime[0]:
            params["structural_elements"] = [tr_5prime[0], gene_3prime[0]]
        elif not tr_5prime[0] and tr_3prime[0]:
            params["structural_elements"] = [gene_5prime[0], tr_3prime[0]]
        else:
            params["structural_elements"] = [tr_5prime[0], tr_3prime[0]]

        if ce:
            params["causative_event"] = ce
        if rf:
            params["r_frame_preserved"] = rf
        return AssayedFusion(**params), None

    def _get_causative_event(
        self, chrom1: str, chrom2: str, descr: Optional[str] = None
    ) -> CausativeEvent:
        """Infer Causative Event. Currently restricted to rearrangements
        :param chrom1: The chromosome for the 5' partner
        :param chrom2: The chromosome for the 3' partner
        :param annots: An annotation describing the fusion event. This input is supplied to the event_description CausativeEvent attribute.
        :return: A CausativeEvent object if construction is successful
        """
        if chrom1 != chrom2 and descr:
            return CausativeEvent(
                event_type=EventType("rearrangement"), event_description=descr
            )
        if chrom1 != chrom2:
            return CausativeEvent(event_type=EventType("rearrangement"))
        return None

    def _get_gene_element_unnormalized(self, symbol: str, caller: str) -> GeneElement:
        """Return GeneElement when gene symbol cannot be normalized
        :param symbol: A gene symbol for a fusion partner
        :return: A GeneElement object
        """
        return (
            GeneElement(
                type="GeneElement",
                gene_descriptor=GeneDescriptor(
                    id=f"gene:{symbol}",
                    type="GeneDescriptor",
                    label=symbol,
                    gene_id=f"{caller}:N/A",
                ),
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
            return (
                ge if ge[0] else self._get_gene_element_unnormalized(genelist, caller)
            )

        genes = genelist.split(",")
        dists = []
        for gene in genes:
            start, end = gene.rfind("("), gene.rfind(")")
            dists.append(int(gene[start + 1 : end]))
        gene = (
            genes[0].split("(")[0] if dists[0] <= dists[1] else genes[1].split("(")[0]
        )
        ge = self.fusor.gene_element(gene=gene)
        return ge if ge[0] else self._get_gene_element_unnormalized(gene, caller)

    async def from_jaffa(self, jaffa_row: pd.DataFrame) -> AssayedFusion:
        """Parse JAFFA fusion output to create AssayedFusion object
        :param jaffa_row: A row of JAFFA output
        :return: An AssayedFusion object, if construction is successful
        """
        genes = jaffa_row["fusion genes"].split(":")
        gene_5prime = self._get_gene_element(genes[0], "jaffa")
        gene_3prime = self._get_gene_element(genes[1], "jaffa")

        # Reformat chromosome, strand orientation
        chrom1 = jaffa_row.chrom1.strip("chr")
        chrom2 = jaffa_row.chrom2.strip("chr")
        strand1 = Strand.POSITIVE if jaffa_row.strand1 == "+" else Strand.NEGATIVE
        strand2 = Strand.POSITIVE if jaffa_row.strand2 == "+" else Strand.NEGATIVE

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=chrom1,
            strand=strand1,
            end=int(jaffa_row.base1),
            gene=genes[0],
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=chrom2,
            strand=strand2,
            start=int(jaffa_row.base2),
            gene=genes[1],
            get_nearest_transcript_junction=True,
        )

        if jaffa_row.rearrangement == "TRUE":
            ce = CausativeEvent(
                event_type=EventType("rearrangement"),
                event_description=jaffa_row["classification"],
            )
        else:
            ce = None

        rf = bool(jaffa_row.inframe == "TRUE")
        return self._format_fusion(
            gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce, rf
        )

    async def from_star_fusion(self, sf_row: pd.DataFrame) -> AssayedFusion:
        """Parse STAR-Fusion output to create AssayedFusion object
        :param sf_row: A row of STAR-Fusion output
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = sf_row["LeftGene"].split("^")[0]
        gene2 = sf_row["RightGene"].split("^")[0]
        gene_5prime = self._get_gene_element(gene1, "star_fusion")
        gene_3prime = self._get_gene_element(gene2, "star_fusion")

        five_prime = sf_row["LeftBreakpoint"].split(":")
        chrom1 = five_prime[0].strip("chr")
        pos1 = int(five_prime[1])
        strand1 = Strand.POSITIVE if five_prime[2] == "+" else Strand.NEGATIVE

        three_prime = sf_row["RightBreakpoint"].split(":")
        chrom2 = three_prime[0].strip("chr")
        pos2 = int(three_prime[1])
        strand2 = Strand.POSITIVE if three_prime[2] == "+" else Strand.NEGATIVE

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=chrom1,
            strand=strand1,
            end=pos1,
            gene=gene1,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=chrom2,
            strand=strand2,
            start=pos2,
            gene=gene2,
            get_nearest_transcript_junction=True,
        )

        ce = self._get_causative_event(chrom1, chrom2, sf_row["annots"])
        return self._format_fusion(gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce)

    async def from_fusion_catcher(self, fc_row: pd.DataFrame) -> AssayedFusion:
        """Parse FusionCatcher output to create AssayedFusion object
        :param fc_row: A row of FusionCatcher output
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = fc_row["Gene_1_symbol(5end_fusion_partner)"]
        gene2 = fc_row["Gene_2_symbol(3end_fusion_partner)"]
        gene_5prime = self._get_gene_element(gene1, "fusion_catcher")
        gene_3prime = self._get_gene_element(gene2, "fusion_catcher")

        five_prime = fc_row["Fusion_point_for_gene_1(5end_fusion_partner)"].split(":")
        chrom1 = five_prime[0]
        pos1 = int(five_prime[1])
        strand1 = Strand.POSITIVE if five_prime[2] == "+" else Strand.NEGATIVE

        three_prime = fc_row["Fusion_point_for_gene_2(3end_fusion_partner)"].split(":")
        chrom2 = three_prime[0]
        pos2 = int(three_prime[1])
        strand2 = Strand.POSITIVE if three_prime[2] == "+" else Strand.NEGATIVE

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=chrom1,
            strand=strand1,
            end=pos1,
            gene=gene1,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=chrom2,
            strand=strand2,
            start=pos2,
            gene=gene2,
            get_nearest_transcript_junction=True,
        )

        ce = self._get_causative_event(chrom1, chrom2, fc_row["Predicted_effect"])
        return self._format_fusion(gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce)

    async def from_fusion_map(self, fmap_row: pd.DataFrame) -> AssayedFusion:
        """Parse FusionMap output to create FUSOR AssayedFusion object
        :param fmap_row: A row of FusionMap output
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = fmap_row["KnownGene1"]
        gene2 = fmap_row["KnownGene2"]
        gene_5prime = self._get_gene_element(gene1, "fusion_map")
        gene_3prime = self._get_gene_element(gene2, "fusion_map")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=fmap_row["Chromosome1"],
            strand=Strand.POSITIVE if fmap_row["Strand"][0] == "+" else Strand.NEGATIVE,
            end=int(fmap_row["Position1"]),
            gene=gene1,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=fmap_row["Chromosome2"],
            strand=Strand.POSITIVE if fmap_row["Strand"][1] == "+" else Strand.NEGATIVE,
            start=int(fmap_row["Position2"]),
            gene=gene2,
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

    async def from_arriba(self, arriba_row: pd.DataFrame) -> AssayedFusion:
        """Parse Arriba output to create AssayedFusion object
        :param arriba_row: A row of Arriba output
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = arriba_row["#gene1"]
        gene2 = arriba_row["gene2"]

        # Arriba reports two gene symbols if a breakpoint occurs in an intergenic
        # region. We select the gene symbol with the smallest distance from the
        # breakpoint.
        gene_5prime = self._get_gene_element(gene1, "arriba")
        gene_3prime = self._get_gene_element(gene2, "arriba")

        breakpoint1 = arriba_row["breakpoint1"].split(":")
        breakpoint2 = arriba_row["breakpoint2"].split(":")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=breakpoint1[0],
            strand=Strand.POSITIVE
            if arriba_row["strand1(gene/fusion)"][2] == "+"
            else Strand.NEGATIVE,
            end=int(breakpoint1[1]),
            gene=gene1,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=breakpoint2[0],
            strand=Strand.POSITIVE
            if arriba_row["strand2(gene/fusion)"][2] == "+"
            else Strand.NEGATIVE,
            start=int(breakpoint2[1]),
            gene=gene2,
            get_nearest_transcript_junction=True,
        )

        ce = (
            CausativeEvent(
                event_type=EventType("read-through"),
                event_description=arriba_row["confidence"],
            )
            if "read_through" in arriba_row["type"]
            else CausativeEvent(
                event_type=EventType("rearrangement"),
                event_description=arriba_row["confidence"],
            )
        )
        rf = bool(arriba_row["reading_frame"] == "in-frame")
        return self._format_fusion(
            gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce, rf
        )

    async def from_cicero(self, cicero_row: pd.DataFrame) -> AssayedFusion:
        """Parse CICERO output to create AssayedFusion object
        :param cicero_row: A row of CICERO output
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = cicero_row["geneA"]
        gene2 = cicero_row["geneB"]
        gene_5prime = self._get_gene_element(gene1, "cicero")
        gene_3prime = self._get_gene_element(gene2, "cicero")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=cicero_row["chrA"],
            strand=Strand.POSITIVE if cicero_row["ortA"] == "+" else Strand.NEGATIVE,
            end=int(cicero_row["posA"]),
            gene=gene1,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=cicero_row["chrB"],
            strand=Strand.POSITIVE if cicero_row["ortB"] == "+" else Strand.NEGATIVE,
            start=int(cicero_row["posB"]),
            gene=gene2,
            get_nearest_transcript_junction=True,
        )

        if cicero_row["type"] == "read_through":
            ce = CausativeEvent(
                event_type=EventType("read-through"),
                event_description=cicero_row["type"],
            )
        else:
            if cicero_row["type"] != "Internal_dup":
                ce = CausativeEvent(
                    event_type=EventType("rearrangement"),
                    event_description=cicero_row["type"],
                )
        rf = bool(int(cicero_row["frame"]) == 1)
        return self._format_fusion(
            gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce, rf
        )

    async def from_mapsplice(self, mapsplice_row: pd.DataFrame) -> AssayedFusion:
        """Parse MapSplice output to create AssayedFusion object
        :param mapsplice_row: A row of MapSplice output
        :retun: An AssayedFusion object, if construction is successful
        """
        gene1 = mapsplice_row[60].strip(",")
        gene2 = mapsplice_row[61].strip(",")
        gene_5prime = self._get_gene_element(gene1, "mapsplice")
        gene_3prime = self._get_gene_element(gene1, "mapsplice")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=mapsplice_row[0].split("~")[0],
            strand=Strand.POSITIVE if mapsplice_row[5][0] == "+" else Strand.NEGATIVE,
            end=int(mapsplice_row[1]),
            gene=gene1,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=mapsplice_row[0].split("~")[1],
            strand=Strand.POSITIVE if mapsplice_row[5][1] == "+" else Strand.NEGATIVE,
            start=int(mapsplice_row[2]),
            gene=gene2,
            get_nearest_transcript_junction=True,
        )

        ce = self._get_causative_event(
            mapsplice_row[0].split("~")[0], mapsplice_row[0].split("~")[1]
        )
        return self._format_fusion(gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce)

    async def from_enfusion(self, enfusion_row: pd.DataFrame) -> AssayedFusion:
        """Parse EnFusion output to create AssayedFusion object
        :param enfusion_row: A row of EnFusion output
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = enfusion_row["Gene1"]
        gene2 = enfusion_row["Gene2"]
        gene_5prime = self._get_gene_element(gene1, "enfusion")
        gene_3prime = self._get_gene_element(gene2, "enfusion")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=enfusion_row["Chr1"],
            strand=Strand.POSITIVE
            if enfusion_row["Strand1"] == "+"
            else Strand.NEGATIVE,
            end=int(enfusion_row["Break1"]),
            gene=gene1,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=enfusion_row["Chr2"],
            strand=Strand.POSITIVE
            if enfusion_row["Strand2"] == "+"
            else Strand.NEGATIVE,
            start=int(enfusion_row["Break2"]),
            gene=gene2,
            get_nearest_transcript_junction=True,
        )

        ce = self._get_causative_event(enfusion_row["Chr1"], enfusion_row["Chr2"])
        return self._format_fusion(gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce)

    async def from_genie(self, genie_row: pd.DataFrame) -> AssayedFusion:
        """Parse GENIE output to create AssayedFusion object
        :param genie_row: A row of EnFusion output
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = genie_row["Site1_Hugo_Symbol"]
        gene2 = genie_row["Site2_Hugo_Symbol"]
        gene_5prime = self._get_gene_element(gene1, "genie")
        gene_3prime = self._get_gene_element(gene1, "genie")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=genie_row["Site1_Chromosome"],
            strand=Strand.POSITIVE
            if "+" in genie_row["Site1_Description"]
            else Strand.NEGATIVE,
            end=int(genie_row["Site1_Position"]),
            gene=gene1,
            get_nearest_transcript_junction=True,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=genie_row["Site2_Chromosome"],
            strand=Strand.POSITIVE
            if "+" in genie_row["Site2_Description"]
            else Strand.NEGATIVE,
            start=int(genie_row["Site2_Position"]),
            gene=gene2,
            get_nearest_transcript_junction=True,
        )

        ce = self._get_causative_event(
            genie_row["Site1_Chromosome"],
            genie_row["Site2_Chromosome"],
            genie_row["Annotation"],
        )
        rf = bool(genie_row["Site2_Effect_on_Frame"] == "In_frame")
        return self._format_fusion(
            gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce, rf
        )
