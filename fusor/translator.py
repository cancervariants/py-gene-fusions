"""Module for translating output from fusion detection methods to FUSOR AssayedFusion
objects
"""

from typing import Optional

import pandas as pd

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
    ) -> AssayedFusion:
        """Format classes to create AssayedFusion objects
        :param gene_5prime: 5'prime GeneElement
        :param gene_3prime: 3'prime GeneElement
        :param Optional[tr_5prime]: 5'prime TranscriptSegmentElement
        :param Optional[tr_3prime]: 3'prime TranscriptSegmentElement
        :param Optional[ce]: CausativeEvent
        :return AssayedFusion object
        """
        if ce is not None:
            if tr_5prime[0] is None and tr_3prime[0] is None:
                return (
                    self.fusor.assayed_fusion(
                        structural_elements=[gene_5prime[0], gene_3prime[0]],
                        causative_event=ce,
                    ),
                    None,
                )
            elif tr_5prime[0] is not None and tr_3prime[0] is None:
                return (
                    self.fusor.assayed_fusion(
                        structural_elements=[tr_5prime[0], gene_3prime[0]],
                        causative_event=ce,
                    ),
                    None,
                )
            elif tr_5prime[0] is None and tr_3prime[0] is not None:
                return (
                    self.fusor.assayed_fusion(
                        structural_elements=[gene_5prime[0], tr_3prime[0]],
                        causative_event=ce,
                    ),
                    None,
                )
            else:
                return (
                    self.fusor.assayed_fusion(
                        structural_elements=[tr_5prime[0], tr_3prime[0]],
                        causative_event=ce,
                    ),
                    None,
                )
        else:
            if tr_5prime[0] is None and tr_3prime[0] is None:
                return (
                    self.fusor.assayed_fusion(
                        structural_elements=[gene_5prime[0], gene_3prime[0]]
                    ),
                    None,
                )
            elif tr_5prime[0] is not None and tr_3prime[0] is None:
                return (
                    self.fusor.assayed_fusion(
                        structural_elements=[tr_5prime[0], gene_3prime[0]]
                    ),
                    None,
                )
            elif tr_5prime[0] is None and tr_3prime[0] is not None:
                return (
                    self.fusor.assayed_fusion(
                        structural_elements=[gene_5prime[0], tr_3prime[0]]
                    ),
                    None,
                )
            else:
                return (
                    self.fusor.assayed_fusion(
                        structural_elements=[tr_5prime[0], tr_3prime[0]]
                    ),
                    None,
                )

    def _get_causative_event(
        self, chrom1: str, chrom2: str, descr: Optional[str] = None
    ) -> CausativeEvent:
        """Infer Causative Event. Currently restricted to rearrangements
        :param chrom1: The chromosome for the 5' partner
        :param chrom2: The chromosome for the 3' partner
        :param annots: An annotation describing the fusion event. This input is supplied to the event_description CausativeEvent attribute.
        :return: A CausativeEvent object if construction is successful
        """
        if chrom1 != chrom2 and descr is not None:
            return CausativeEvent(
                event_type=EventType("rearrangement"), event_description=descr
            )
        elif chrom1 != chrom2:
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
                    id=f"fusor.gene:{symbol}",
                    type="GeneDescriptor",
                    label=symbol,
                    gene_id=f"{caller}:N/A",
                ),
            ),
            None,
        )

    async def from_jaffa(self, jaffa_row: pd.DataFrame) -> AssayedFusion:
        """Parse JAFFA fusion output to create AssayedFusion object
        :param jaffa_row: A row of JAFFA output
        :return: An AssayedFusion object, if construction is successful
        """
        genes = jaffa_row["fusion genes"].split(":")
        gene_5prime = self.fusor.gene_element(gene=genes[0])
        gene_3prime = self.fusor.gene_element(gene=genes[1])
        if gene_5prime[0] is None:
            gene_5prime = self._get_gene_element_unnormalized(genes[0], "jaffa_5prime")
        if gene_3prime[0] is None:
            gene_3prime = self._get_gene_element_unnormalized(genes[1], "jaffa_3prime")

        # Reformat chromosome, strand orientation
        chrom1 = jaffa_row.chrom1.strip("chr")
        chrom2 = jaffa_row.chrom2.strip("chr")
        strand1 = 1 if jaffa_row.strand1 == "+" else -1
        strand2 = 1 if jaffa_row.strand1 == "+" else -1

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=chrom1,
            strand=strand1,
            end=int(jaffa_row.base1),
            gene=genes[0],
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=chrom2,
            strand=strand2,
            start=int(jaffa_row.base2),
            gene=genes[1],
        )

        if jaffa_row.rearrangement:
            ce = CausativeEvent(
                event_type=EventType("rearrangement"),
                event_description=jaffa_row["classification"],
            )
        else:
            ce = None

        return self._format_fusion(gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce)

    async def from_star_fusion(self, sf_row: pd.DataFrame) -> AssayedFusion:
        """Parse STAR-Fusion output to create AssayedFusion object
        :param sf_row: A row of STAR-Fusion output
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = sf_row["LeftGene"].split("^")[0]
        gene2 = sf_row["RightGene"].split("^")[0]
        gene_5prime = self.fusor.gene_element(gene=gene1)
        gene_3prime = self.fusor.gene_element(gene=gene2)
        if gene_5prime[0] is None:
            gene_5prime = self._get_gene_element_unnormalized(
                gene1, "star_fusion_5prime"
            )
        if gene_3prime[0] is None:
            gene_3prime = self._get_gene_element_unnormalized(
                gene2, "star_fusion_3prime"
            )

        five_prime = sf_row["LeftBreakpoint"].split(":")
        chrom1 = five_prime[0].strip("chr")
        pos1 = int(five_prime[1])
        strand1 = 1 if five_prime[2] == "+" else -1

        three_prime = sf_row["RightBreakpoint"].split(":")
        chrom2 = three_prime[0].strip("chr")
        pos2 = int(three_prime[1])
        strand2 = 1 if three_prime[2] == "+" else -1

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=chrom1,
            strand=strand1,
            end=pos1,
            gene=gene1,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=chrom2,
            strand=strand2,
            start=pos2,
            gene=gene2,
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
        gene_5prime = self.fusor.gene_element(gene=gene1)
        gene_3prime = self.fusor.gene_element(gene=gene2)
        if gene_5prime[0] is None:
            gene_5prime = self._get_gene_element_unnormalized(
                gene1, "fusion_catcher_5prime"
            )
        if gene_3prime[0] is None:
            gene_3prime = self._get_gene_element_unnormalized(
                gene2, "fusion_catcher_3prime"
            )

        five_prime = fc_row["Fusion_point_for_gene_1(5end_fusion_partner)"].split(":")
        chrom1 = five_prime[0]
        pos1 = int(five_prime[1])
        strand1 = 1 if five_prime[2] == "+" else -1

        three_prime = fc_row["Fusion_point_for_gene_2(3end_fusion_partner)"].split(":")
        chrom2 = three_prime[0]
        pos2 = int(three_prime[1])
        strand2 = 1 if three_prime[2] == "+" else -1

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=chrom1,
            strand=strand1,
            end=pos1,
            gene=gene1,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=chrom2,
            strand=strand2,
            start=pos2,
            gene=gene2,
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
        gene_5prime = self.fusor.gene_element(gene=gene1)
        gene_3prime = self.fusor.gene_element(gene=gene2)
        if gene_5prime[0] is None:
            gene_5prime = self._get_gene_element_unnormalized(
                gene1, "fusion_map_5prime"
            )
        if gene_3prime[0] is None:
            gene_3prime = self._get_gene_element_unnormalized(
                gene2, "fusion_map_3prime"
            )

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=fmap_row["Chromosome1"],
            strand=1 if fmap_row["Strand"][0] == "+" else -1,
            end=int(fmap_row["Position1"]),
            gene=gene1,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=fmap_row["Chromosome2"],
            strand=1 if fmap_row["Strand"][1] == "+" else -1,
            start=int(fmap_row["Position2"]),
            gene=gene2,
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
        return self._format_fusion(gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce)

    async def from_arriba(self, arriba_row: pd.DataFrame) -> AssayedFusion:
        """Parse Arriba output to create AssayedFusion object
        :param arriba_row: A row of Arriba output
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = arriba_row["gene1"]
        gene2 = arriba_row["gene2"]
        gene_5prime = self.fusor.gene_element(gene=gene1)
        gene_3prime = self.fusor.gene_element(gene=gene2)
        if gene_5prime[0] is None:
            gene_5prime = self._get_gene_element_unnormalized(gene1, "arriba_5prime")
        if gene_3prime[0] is None:
            gene_3prime = self._get_gene_element_unnormalized(gene2, "arriba_3prime")

        breakpoint1 = arriba_row["breakpoint1"].split(":")
        breakpoint2 = arriba_row["breakpoint2"].split(":")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=breakpoint1[0],
            strand=1 if arriba_row["strand1(gene/fusion)"][2] == "+" else -1,
            end=int(breakpoint1[1]),
            gene=gene1,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=breakpoint2[0],
            strand=1 if arriba_row["strand2(gene/fusion)"][2] == "+" else -1,
            start=int(breakpoint2[1]),
            gene=gene2,
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
        return self._format_fusion(gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce)

    async def from_cicero(self, cicero_row: pd.DataFrame) -> AssayedFusion:
        """Parse CICERO output to create AssayedFusion object
        :param cicero_row: A row of CICERO output
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = cicero_row["geneA"]
        gene2 = cicero_row["geneB"]
        gene_5prime = self.fusor.gene_element(gene=gene1)
        gene_3prime = self.fusor.gene_element(gene=gene2)
        if gene_5prime[0] is None:
            gene_5prime = self._get_gene_element_unnormalized(gene1, "cicero_5prime")
        if gene_3prime[0] is None:
            gene_3prime = self._get_gene_element_unnormalized(gene2, "cicero_3prime")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=cicero_row["chrA"],
            strand=1 if cicero_row["ortA"] == "+" else -1,
            end=int(cicero_row["posA"]),
            gene=gene1,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=cicero_row["chrB"],
            strand=1 if cicero_row["ortB"] == "+" else -1,
            start=int(cicero_row["posB"]),
            gene=gene2,
        )

        if cicero_row["type"] == "read_through":
            ce = CausativeEvent(
                event_type=EventType("read-through"),
                event_description=cicero_row["rating"],
            )
        else:
            if cicero_row["type"] != "Internal_dup":
                ce = CausativeEvent(
                    event_type=EventType("rearrangement"),
                    event_description=cicero_row["rating"],
                )
        return self._format_fusion(gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce)

    async def from_mapsplice(self, mapsplice_row: pd.DataFrame) -> AssayedFusion:
        """Parse MapSplice output to create AssayedFusion object
        :param mapsplice_row: A row of MapSplice output
        :retun: An AssayedFusion object, if construction is successful
        """
        gene1 = mapsplice_row[60].strip(",")
        gene2 = mapsplice_row[61].strip(",")
        gene_5prime = self.fusor.gene_element(gene=gene1)
        gene_3prime = self.fusor.gene_element(gene=gene2)
        if gene_5prime[0] is None:
            gene_5prime = self._get_gene_element_unnormalized(gene1, "mapslice_5prime")
        if gene_3prime[0] is None:
            gene_3prime = self._get_gene_element_unnormalized(gene2, "mapslice_3prime")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=mapsplice_row[0].split("~")[0],
            strand=1 if mapsplice_row[5][0] == "+" else -1,
            end=int(mapsplice_row[1]),
            gene=gene1,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=mapsplice_row[0].split("~")[1],
            strand=1 if mapsplice_row[5][1] == "+" else -1,
            start=int(mapsplice_row[2]),
            gene=gene2,
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
        gene_5prime = self.fusor.gene_element(gene=gene1)
        gene_3prime = self.fusor.gene_element(gene=gene2)
        if gene_5prime[0] is None:
            gene_5prime = self._get_gene_element_unnormalized(gene1, "enfusion_5prime")
        if gene_3prime[0] is None:
            gene_3prime = self._get_gene_element_unnormalized(gene2, "enufusion_3prime")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=enfusion_row["Chr1"],
            strand=1 if enfusion_row["Strand1"] == "+" else -1,
            end=int(enfusion_row["Break1"]),
            gene=gene1,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=enfusion_row["Chr2"],
            strand=1 if enfusion_row["Strand2"] == "+" else -1,
            start=int(enfusion_row["Break2"]),
            gene=gene2,
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
        gene_5prime = self.fusor.gene_element(gene=gene1)
        gene_3prime = self.fusor.gene_element(gene=gene2)
        if gene_5prime[0] is None:
            gene_5prime = self._get_gene_element_unnormalized(gene1, "genie_5prime")
        if gene_3prime[0] is None:
            gene_3prime = self._get_gene_element_unnormalized(gene2, "genie_3prime")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=genie_row["Site1_Chromosome"],
            strand=1 if "+" in genie_row["Site1_Description"] else -1,
            end=int(genie_row["Site1_Position"]),
            gene=gene1,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            chromosome=genie_row["Site2_Chromosome"],
            strand=1 if "+" in genie_row["Site2_Description"] else -1,
            start=int(genie_row["Site2_Position"]),
            gene=gene2,
        )

        ce = self._get_causative_event(
            genie_row["Site1_Chromosome"],
            genie_row["Site2_Chromosome"],
            genie_row["Annotation"],
        )
        return self._format_fusion(gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce)
