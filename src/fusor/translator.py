"""Module for translating output from fusion detection methods to FUSOR AssayedFusion
objects
"""

import logging

import polars as pl
from cool_seq_tool.schemas import Assembly, CoordinateType
from pydantic import BaseModel

from fusor.fusion_caller_models import (
    JAFFA,
    Arriba,
    Caller,
    Cicero,
    EnFusion,
    FusionCatcher,
    Genie,
    STARFusion,
)
from fusor.fusor import FUSOR
from fusor.models import (
    AnchoredReads,
    Assay,
    AssayedFusion,
    BreakpointCoverage,
    CausativeEvent,
    ContigSequence,
    EventType,
    GeneElement,
    ReadData,
    SpanningReads,
    SplitReads,
    TranscriptSegmentElement,
    UnknownGeneElement,
)

_logger = logging.getLogger(__name__)


class GeneFusionPartners(BaseModel):
    """Class for defining gene fusion partners"""

    gene_5prime_element: GeneElement | UnknownGeneElement
    gene_5prime: str | None = None
    gene_3prime_element: GeneElement | UnknownGeneElement
    gene_3prime: str | None = None


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
        gene_5prime: GeneElement | UnknownGeneElement,
        gene_3prime: GeneElement | UnknownGeneElement,
        tr_5prime: TranscriptSegmentElement | None = None,
        tr_3prime: TranscriptSegmentElement | None = None,
        ce: CausativeEvent | None = None,
        rf: bool | None = None,
        assay: Assay | None = None,
        contig: ContigSequence | None = None,
        reads: ReadData | None = None,
    ) -> AssayedFusion:
        """Format classes to create AssayedFusion objects

        :param gene_5prime: 5'prime GeneElement
        :param gene_3prime: 3'prime GeneElement
        :param tr_5prime: 5'prime TranscriptSegmentElement
        :param tr_3prime: 3'prime TranscriptSegmentElement
        :param ce: CausativeEvent
        :param rf: A boolean indicating if the reading frame is preserved
        :param assay: Assay
        :param contig: The contig sequence
        :param reads: The read data
        :return AssayedFusion object
        """
        params = {
            "causativeEvent": ce,
            "readingFramePreserved": rf,
            "assay": assay,
            "contig": contig,
            "readData": reads,
        }
        if not tr_5prime[0] and not tr_3prime[0]:
            params["structure"] = [gene_5prime, gene_3prime]
        elif tr_5prime[0] and not tr_3prime[0]:
            params["structure"] = [tr_5prime[0], gene_3prime]
        elif not tr_5prime[0] and tr_3prime[0]:
            params["structure"] = [gene_5prime, tr_3prime[0]]
        else:
            params["structure"] = [tr_5prime[0], tr_3prime[0]]
        return AssayedFusion(**params)

    def _get_causative_event(
        self, chrom1: str, chrom2: str, descr: str | None = None
    ) -> CausativeEvent | None:
        """Infer Causative Event. Currently restricted to rearrangements

        :param chrom1: The chromosome for the 5' partner
        :param chrom2: The chromosome for the 3' partner
        :param descr: An annotation describing the fusion event. This input is supplied to the eventDescription CausativeEvent attribute.
        :return: A CausativeEvent object if construction is successful
        """
        if descr and "rearrangement" in descr:
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
        return GeneElement(
            gene={
                "id": f"gene:{symbol}",
                "label": symbol,
                "type": "Gene",
            },
        )

    def _get_gene_element(self, genes: str, caller: Caller) -> GeneElement:
        """Return a GeneElement given an individual/list of gene symbols and a
        fusion detection algorithm

        :param genes: A gene symbol or list of gene symbols, separated by columns
        :param caller: The examined fusion detection algorithm
        :return A GeneElement object
        """
        if "," not in genes or caller != caller.ARRIBA:
            ge = self.fusor.gene_element(gene=genes)
            return ge[0] if ge[0] else self._get_gene_element_unnormalized(genes)

        genes = genes.split(",")
        dists = []
        for gene in genes:
            start, end = gene.rfind("("), gene.rfind(")")
            dists.append(int(gene[start + 1 : end]))
        gene = (
            genes[0].split("(")[0] if dists[0] <= dists[1] else genes[1].split("(")[0]
        )
        ge = self.fusor.gene_element(gene=gene)
        return ge[0] if ge[0] else self._get_gene_element_unnormalized(gene)

    def _are_fusion_partners_different(
        self,
        gene_5prime: str | UnknownGeneElement,
        gene_3prime: str | UnknownGeneElement,
    ) -> bool:
        """Check if the normalized gene symbols for the two fusion partners
        are different. If not, this event is not a fusion

        :param gene_5prime: The 5' gene partner
        :param gene_3prime: The 3' gene partner
        :return ``True`` if the gene symbols are different, ``False`` if not
        """
        if gene_5prime != gene_3prime:
            return True
        _logger.error(
            "The supplied fusion is not valid as the two fusion partners are the same"
        )
        return False

    def _get_genomic_ac(self, chrom: str, build: Assembly) -> str:
        """Return a RefSeq genomic accession given a chromosome and a reference build

        :param chrom: A chromosome number
        :param build: The assembly, either GRCh37 or GRCh38
        :return: The corresponding refseq genomic accession
        :raise ValueError: if unable to retrieve genomic accession
        """
        sr = self.fusor.cool_seq_tool.seqrepo_access
        alias_list, errors = sr.translate_identifier(
            f"{build}:{chrom}", target_namespaces="refseq"
        )
        if errors:
            statement = f"Genomic accession for {chrom} could not be retrieved"
            _logger.error(statement)
            raise ValueError
        return alias_list[0].split(":")[1]

    def _process_gene_symbols(
        self, gene_5prime: str, gene_3prime: str, caller: Caller
    ) -> dict[GeneElement | UnknownGeneElement]:
        """Process gene symbols to create GeneElements or UnknownGeneElements

        :param gene_5prime: The 5' gene symbol
        :param gene_3prime: The 3' gene symbol
        :param caller: The gene fusion caller
        :return A dictionary of GeneElements or UnknownGeneElements
        """
        gene_5prime_element = gene_3prime_element = None
        if gene_5prime == "NA":
            gene_5prime_element = gene_5prime = UnknownGeneElement()
            gene_5prime = None
        if gene_3prime == "NA":
            gene_3prime_element = gene_3prime = UnknownGeneElement()
            gene_3prime = None
        if not gene_5prime_element:
            gene_5prime_element = self._get_gene_element(gene_5prime, caller)
            gene_5prime = gene_5prime_element.gene.label
        if not gene_3prime_element:
            gene_3prime_element = self._get_gene_element(gene_3prime, caller)
            gene_3prime = gene_3prime_element.gene.label
        params = {
            "gene_5prime_element": gene_5prime_element,
            "gene_5prime": gene_5prime,
            "gene_3prime_element": gene_3prime_element,
            "gene_3prime": gene_3prime,
        }
        return GeneFusionPartners(**params)

    async def from_jaffa(
        self,
        jaffa: JAFFA,
        coordinate_type: CoordinateType,
        rb: Assembly,
    ) -> AssayedFusion | None:
        """Parse JAFFA fusion output to create AssayedFusion object

        :param JAFFA: Output from JAFFA caller
        :param coordinate_type: If the coordinate is inter-residue or residue
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        genes = jaffa.fusion_genes.split(":")
        fusion_partners = self._process_gene_symbols(genes[0], genes[1], Caller.JAFFA)

        if not self._are_fusion_partners_different(
            fusion_partners.gene_5prime, fusion_partners.gene_3prime
        ):
            return None

        if not isinstance(fusion_partners.gene_5prime_element, UnknownGeneElement):
            tr_5prime = await self.fusor.transcript_segment_element(
                tx_to_genomic_coords=False,
                genomic_ac=self._get_genomic_ac(jaffa.chrom1, rb),
                seg_end_genomic=jaffa.base1,
                gene=fusion_partners.gene_5prime,
                coordinate_type=coordinate_type,
                starting_assembly=rb,
            )

        if not isinstance(fusion_partners.gene_3prime_element, UnknownGeneElement):
            tr_3prime = await self.fusor.transcript_segment_element(
                tx_to_genomic_coords=False,
                genomic_ac=self._get_genomic_ac(jaffa.chrom2, rb),
                seg_start_genomic=jaffa.base2,
                gene=fusion_partners.gene_3prime,
                coordinate_type=coordinate_type,
                starting_assembly=rb,
            )

        if jaffa.rearrangement:
            ce = CausativeEvent(
                eventType=EventType("rearrangement"),
                eventDescription=jaffa.classification,
            )
        else:
            ce = None

        read_data = ReadData(
            split=SplitReads(splitReads=jaffa.spanning_reads),
            spanning=SpanningReads(spanningReads=jaffa.spanning_pairs),
        )

        return self._format_fusion(
            fusion_partners.gene_5prime_element,
            fusion_partners.gene_3prime_element,
            tr_5prime
            if isinstance(fusion_partners.gene_5prime_element, GeneElement)
            else [None],
            tr_3prime
            if isinstance(fusion_partners.gene_3prime_element, GeneElement)
            else [None],
            ce,
            jaffa.inframe if isinstance(jaffa.inframe, bool) else None,
            reads=read_data,
        )

    async def from_star_fusion(
        self,
        star_fusion: STARFusion,
        coordinate_type: CoordinateType,
        rb: Assembly,
    ) -> AssayedFusion:
        """Parse STAR-Fusion output to create AssayedFusion object

        :param star_fusion: Output from STAR-Fusion caller
        :param coordinate_type: If the coordinate is inter-residue or residue
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = star_fusion.left_gene.split("^")[0]
        gene2 = star_fusion.right_gene.split("^")[0]

        fusion_partners = self._process_gene_symbols(gene1, gene2, Caller.STAR_FUSION)
        if not self._are_fusion_partners_different(
            fusion_partners.gene_5prime, fusion_partners.gene_3prime
        ):
            return None

        five_prime = star_fusion.left_breakpoint.split(":")
        three_prime = star_fusion.right_breakpoint.split(":")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(five_prime[0], rb),
            seg_end_genomic=int(five_prime[1]),
            gene=fusion_partners.gene_5prime,
            coordinate_type=coordinate_type,
            starting_assembly=rb,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(three_prime[0], rb),
            seg_start_genomic=int(three_prime[1]),
            gene=fusion_partners.gene_3prime,
            coordinate_type=coordinate_type,
            starting_assembly=rb,
        )

        ce = self._get_causative_event(
            five_prime[0], three_prime[0], ",".join(star_fusion.annots)
        )
        read_data = ReadData(
            split=SplitReads(splitReads=star_fusion.junction_read_count),
            spanning=SpanningReads(spanningReads=star_fusion.spanning_frag_count),
        )

        return self._format_fusion(
            fusion_partners.gene_5prime_element,
            fusion_partners.gene_3prime_element,
            tr_5prime
            if isinstance(fusion_partners.gene_5prime_element, GeneElement)
            else [None],
            tr_3prime
            if isinstance(fusion_partners.gene_3prime_element, GeneElement)
            else [None],
            ce,
            reads=read_data,
        )

    async def from_fusion_catcher(
        self,
        fusion_catcher: FusionCatcher,
        coordinate_type: CoordinateType,
        rb: Assembly,
    ) -> AssayedFusion:
        """Parse FusionCatcher output to create AssayedFusion object

        :param fusion_catcher: Output from FusionCatcher caller
        :param coordinate_type: If the coordinate is inter-residue or residue
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene_5prime_element = self._get_gene_element(
            fusion_catcher.five_prime_partner, Caller.FUSION_CATCHER
        )
        gene_3prime_element = self._get_gene_element(
            fusion_catcher.three_prime_partner, Caller.FUSION_CATCHER
        )
        if not self._are_fusion_partners_different(
            gene_5prime_element.gene.label, gene_3prime_element.gene.label
        ):
            return None

        five_prime = fusion_catcher.five_prime_fusion_point.split(":")
        three_prime = fusion_catcher.three_prime_fusion_point.split(":")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(five_prime[0], rb),
            seg_end_genomic=int(five_prime[1]),
            gene=gene_5prime_element.gene.label,
            coordinate_type=coordinate_type,
            starting_assembly=rb,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(three_prime[0], rb),
            seg_start_genomic=int(three_prime[1]),
            gene=gene_3prime_element.gene.label,
            coordinate_type=coordinate_type,
            starting_assembly=rb,
        )

        ce = self._get_causative_event(
            five_prime[0], three_prime[0], fusion_catcher.predicted_effect
        )
        read_data = ReadData(
            split=SplitReads(splitReads=fusion_catcher.spanning_unique_reads),
            spanning=SpanningReads(spanningReads=fusion_catcher.spanning_reads),
        )
        contig = ContigSequence(contig=fusion_catcher.fusion_sequence)

        return self._format_fusion(
            gene_5prime_element,
            gene_3prime_element,
            tr_5prime,
            tr_3prime,
            ce,
            contig=contig,
            reads=read_data,
        )

    async def from_fusion_map(
        self, fmap_row: pl.DataFrame, coordinate_type: CoordinateType, rb: Assembly
    ) -> AssayedFusion:
        """Parse FusionMap output to create FUSOR AssayedFusion object

        :param fmap_row: A row of FusionMap output
        :param rb: The reference build used to call the fusion
        :param coordinate_type: If the coordinate is inter-residue or residue
        :return: An AssayedFusion object, if construction is successful
        """
        gene1 = fmap_row.get_column("KnownGene1").item()
        gene2 = fmap_row.get_column("KnownGene2").item()
        gene_5prime = self._get_gene_element(gene1, "fusion_map").gene.label
        gene_3prime = self._get_gene_element(gene2, "fusion_map").gene.label

        if not self._are_fusion_partners_different(gene_5prime, gene_3prime):
            return None

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(
                fmap_row.get_column("Chromosome1").item(), rb
            ),
            seg_end_genomic=int(fmap_row.get_column("Position1").item()),
            gene=gene_5prime,
            coordinate_type=coordinate_type,
            starting_assembly=rb,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(
                fmap_row.get_column("Chromosome2").item(), rb
            ),
            seg_start_genomic=int(fmap_row.get_column("Position2").item()),
            gene=gene_3prime,
            coordinate_type=coordinate_type,
            starting_assembly=rb,
        )

        # Combine columns to create fusion annotation string"
        descr = (
            fmap_row.get_column("FusionGene").item()
            + ","
            + fmap_row.get_column("SplicePatternClass").item()
            + ","
            + fmap_row.get_column("FrameShiftClass").item()
        )
        ce = self._get_causative_event(
            fmap_row.get_column("Chromosome1").item(),
            fmap_row.get_column("Chromosome2").item(),
            descr,
        )
        rf = bool(fmap_row.get_column("FrameShiftClass").item() == "InFrame")
        return self._format_fusion(
            gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce, rf
        )

    async def from_arriba(
        self,
        arriba: Arriba,
        coordinate_type: CoordinateType,
        rb: Assembly,
    ) -> AssayedFusion:
        """Parse Arriba output to create AssayedFusion object

        :param arriba: Output from Arriba caller
        :param coordinate_type: If the coordinate is inter-residue or residue
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        # Arriba reports two gene symbols if a breakpoint occurs in an intergenic
        # space. We select the gene symbol with the smallest distance from the
        # breakpoint.
        gene_5prime_element = self._get_gene_element(arriba.gene1, "arriba")
        gene_3prime_element = self._get_gene_element(arriba.gene2, "arriba")
        gene_5prime = gene_5prime_element.gene.label
        gene_3prime = gene_3prime_element.gene.label

        if not self._are_fusion_partners_different(gene_5prime, gene_3prime):
            return None

        strand1 = arriba.strand1.split("/")[1]  # Determine strand that is transcribed
        strand2 = arriba.strand2.split("/")[1]  # Determine strand that is transcribed
        if strand1 == "+":
            gene1_seg_start = arriba.direction1 == "upstream"
        else:
            gene1_seg_start = arriba.direction1 == "downstream"
        if strand2 == "+":
            gene2_seg_start = arriba.direction2 == "upstream"
        else:
            gene2_seg_start = arriba.direction2 == "downstream"

        breakpoint1 = arriba.breakpoint1.split(":")
        breakpoint2 = arriba.breakpoint2.split(":")

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(breakpoint1[0], rb),
            seg_start_genomic=int(breakpoint1[1]) if gene1_seg_start else None,
            seg_end_genomic=int(breakpoint1[1]) if not gene1_seg_start else None,
            gene=gene_5prime,
            coverage=BreakpointCoverage(fragmentCoverage=arriba.coverage1),
            reads=AnchoredReads(reads=arriba.split_reads1),
            coordinate_type=coordinate_type,
            starting_assembly=rb,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(breakpoint2[0], rb),
            seg_start_genomic=int(breakpoint2[1]) if gene2_seg_start else None,
            seg_end_genomic=int(breakpoint2[1]) if not gene2_seg_start else None,
            gene=gene_3prime,
            coverage=BreakpointCoverage(fragmentCoverage=arriba.coverage2),
            reads=AnchoredReads(reads=arriba.split_reads2),
            coordinate_type=coordinate_type,
            starting_assembly=rb,
        )

        ce = (
            CausativeEvent(
                eventType=EventType("read-through"),
                eventDescription=arriba.confidence,
            )
            if "read_through" in arriba.event_type
            else CausativeEvent(
                eventType=EventType("rearrangement"),
                eventDescription=arriba.confidence,
            )
        )
        rf = bool(arriba.rf == "in-frame") if arriba.rf != "." else None
        read_data = ReadData(
            spanning=SpanningReads(spanningReads=arriba.discordant_mates)
        )
        contig = ContigSequence(contig=arriba.fusion_transcript)

        return self._format_fusion(
            gene_5prime_element,
            gene_3prime_element,
            tr_5prime,
            tr_3prime,
            ce,
            rf,
            contig=contig,
            reads=read_data,
        )

    async def from_cicero(
        self,
        cicero: Cicero,
        coordinate_type: CoordinateType,
        rb: Assembly,
    ) -> AssayedFusion | str:
        """Parse CICERO output to create AssayedFusion object

        :param cicero: Output from CICERO caller
        :param coordinate_type: If the coordinate is inter-residue or residue
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        # Check if gene symbols have valid formatting. CICERO can output two or more
        # gene symbols for `gene_5prime` or `gene_3prime`, which are separated by a comma. As
        # there is not a precise way to resolve this ambiguity, we do not process
        # these events
        if "," in cicero.gene_5prime or "," in cicero.gene_3prime:
            msg = "Ambiguous gene symbols are reported by CICERO for at least one of the fusion partners"
            _logger.warning(msg)
            return msg

        # Check CICERO annotation regarding the confidence that the called fusion
        # has biological meaning
        if cicero.sv_ort != ">":
            msg = "CICERO annotation indicates that this event does not have confident biological meaning"
            _logger.warning(msg)
            return msg

        gene_5prime_element = self._get_gene_element(cicero.gene_5prime, "cicero")
        gene_3prime_element = self._get_gene_element(cicero.gene_3prime, "cicero")
        gene_5prime = gene_5prime_element.gene.label
        gene_3prime = gene_3prime_element.gene.label

        if not self._are_fusion_partners_different(gene_5prime, gene_3prime):
            return None

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(cicero.chr_5prime, rb),
            seg_end_genomic=cicero.pos_5prime,
            gene=gene_5prime,
            coverage=BreakpointCoverage(fragmentCoverage=cicero.coverage_5prime),
            reads=AnchoredReads(reads=cicero.reads_5prime),
            coordinate_type=coordinate_type,
            starting_assembly=rb,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(cicero.chr_3prime, rb),
            seg_start_genomic=cicero.pos_3prime,
            gene=gene_3prime,
            coverage=BreakpointCoverage(fragmentCoverage=cicero.coverage_3prime),
            reads=AnchoredReads(reads=cicero.reads_3prime),
            coordinate_type=coordinate_type,
            starting_assembly=rb,
        )

        if cicero.event_type == "read_through":
            ce = CausativeEvent(
                eventType=EventType("read-through"),
                eventDescription=cicero.event_type,
            )
        else:
            ce = CausativeEvent(
                eventType=EventType("rearrangement"),
                eventDescription=cicero.event_type,
            )
        contig = ContigSequence(contig=cicero.contig)

        return self._format_fusion(
            gene_5prime_element,
            gene_3prime_element,
            tr_5prime,
            tr_3prime,
            ce,
            contig=contig,
        )

    async def from_mapsplice(
        self, mapsplice_row: pl.DataFrame, coordinate_type: CoordinateType, rb: Assembly
    ) -> AssayedFusion:
        """Parse MapSplice output to create AssayedFusion object

        :param mapsplice_row: A row of MapSplice output
        :param rb: The reference build used to call the fusion
        :param coordinate_type: If the coordinate is inter-residue or residue
        :retun: An AssayedFusion object, if construction is successful
        """
        gene1 = mapsplice_row[60].strip(",")
        gene2 = mapsplice_row[61].strip(",")
        gene_5prime_element = self._get_gene_element(gene1, "mapsplice")
        gene_3prime_element = self._get_gene_element(gene2, "mapsplice")
        gene_5prime = gene_5prime_element.gene.label
        gene_3prime = gene_3prime_element.gene.label

        if not self._are_fusion_partners_different(gene_5prime, gene_3prime):
            return None

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(mapsplice_row[0].split("~")[0], rb),
            seg_end_genomic=int(mapsplice_row[1]),
            gene=gene_5prime,
            coordinate_type=coordinate_type,
            starting_assembly=rb,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(mapsplice_row[0].split("~")[1], rb),
            seg_start_genomic=int(mapsplice_row[2]),
            gene=gene_3prime,
            coordinate_type=coordinate_type,
            starting_assembly=rb,
        )

        ce = self._get_causative_event(
            mapsplice_row[0].split("~")[0], mapsplice_row[0].split("~")[1]
        )
        return self._format_fusion(gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce)

    async def from_enfusion(
        self,
        enfusion: EnFusion,
        coordinate_type: CoordinateType,
        rb: Assembly,
    ) -> AssayedFusion:
        """Parse EnFusion output to create AssayedFusion object

        :param enfusion: Output from EnFusion caller
        :param coordinate_type: If the coordinate is inter-residue or residue
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene_5prime_element = self._get_gene_element(enfusion.gene_5prime, "enfusion")
        gene_3prime_element = self._get_gene_element(enfusion.gene_3prime, "enfusion")
        gene_5prime = gene_5prime_element.gene.label
        gene_3prime = gene_3prime_element.gene.label

        if not self._are_fusion_partners_different(gene_5prime, gene_3prime):
            return None

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(enfusion.chr_5prime, rb),
            seg_end_genomic=enfusion.break_5prime,
            gene=gene_5prime,
            coordinate_type=coordinate_type,
            starting_assembly=rb,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(enfusion.chr_3prime, rb),
            seg_start_genomic=enfusion.break_3prime,
            gene=gene_3prime,
            coordinate_type=coordinate_type,
            starting_assembly=rb,
        )

        ce = self._get_causative_event(
            enfusion.chr_5prime,
            enfusion.chr_3prime,
        )
        return self._format_fusion(
            gene_5prime_element, gene_3prime_element, tr_5prime, tr_3prime, ce
        )

    async def from_genie(
        self,
        genie: Genie,
        coordinate_type: CoordinateType,
        rb: Assembly,
    ) -> AssayedFusion:
        """Parse GENIE output to create AssayedFusion object

        :param genie: Output from GENIE dataset
        :param coordinate_type: If the coordinate is inter-residue or residue
        :param rb: The reference build used to call the fusion
        :return: An AssayedFusion object, if construction is successful
        """
        gene_5prime_element = self._get_gene_element(genie.site1_hugo, "genie")
        gene_3prime_element = self._get_gene_element(genie.site2_hugo, "genie")
        gene_5prime = gene_5prime_element.gene.label
        gene_3prime = gene_3prime_element.gene.label

        if not self._are_fusion_partners_different(gene_5prime, gene_3prime):
            return None

        tr_5prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(genie.site1_chrom, rb),
            seg_end_genomic=genie.site1_pos,
            gene=gene_5prime,
            coordinate_type=coordinate_type,
            starting_assembly=rb,
        )

        tr_3prime = await self.fusor.transcript_segment_element(
            tx_to_genomic_coords=False,
            genomic_ac=self._get_genomic_ac(genie.site2_chrom, rb),
            seg_start_genomic=genie.site2_pos,
            gene=gene_3prime,
            coordinate_type=coordinate_type,
            starting_assembly=rb,
        )

        ce = self._get_causative_event(
            genie.site1_chrom,
            genie.site2_chrom,
            genie.annot,
        )
        rf = bool(genie.reading_frame == "in frame")
        return self._format_fusion(
            gene_5prime, gene_3prime, tr_5prime, tr_3prime, ce, rf
        )
