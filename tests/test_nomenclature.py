"""Test nomenclature generation."""
import json

import pytest

from fusor.models import CategoricalFusion, AssayedFusion
from tests.conftest import EXAMPLES_DIR


@pytest.fixture(scope="module")
def igh_myc():
    """Provide IGH-MYC fusion fixture"""
    return CategoricalFusion(**{
        "structural_components": [
            {
                "type": "GeneComponent",
                "gene_descriptor": {
                    "type": "GeneDescriptor",
                    "gene_id": "hgnc:5477",
                    "label": "IGH",
                    "id": "gene:IGH"
                }
            },
            {
                "type": "GeneComponent",
                "gene_descriptor": {
                    "type": "GeneDescriptor",
                    "gene_id": "hgnc:7553",
                    "label": "MYC",
                    "id": "gene:MYC"
                }
            }
        ],
        "regulatory_elements": [
            {
                "type": "RegulatoryElement",
                "element_type": "promoter",
                "associated_gene": {
                    "type": "GeneDescriptor",
                    "gene_id": "hgnc:5477",
                    "label": "IGH",
                    "id": "gene:IGH"
                }
            }
        ]
    })


def test_generate_nomenclature(fusor, fusion, fusion_example,
                               exhaustive_example, igh_myc):
    """Test that nomenclature generation is correct."""
    fixture_nomenclature = "reg_promoter@BRAF(hgnc:1097)::refseq:NM_152263.3(TPM3):e.1_8::ALK(hgnc:427)::ACGT::refseq:NC_000012.12(chr 12):g.44908821_44908822(+)::v"  # noqa: E501
    nm = fusor.generate_nomenclature(CategoricalFusion(**fusion))
    assert nm == fixture_nomenclature

    fixture_nomenclature = "reg_promoter@BRAF(hgnc:1097)::refseq:NM_152263.3(TPM3):e.1_8::ALK(hgnc:427)::ACGT::refseq:NC_000023.11(chr X):g.44908820_44908822(+)::v"  # noqa: E501
    nm = fusor.generate_nomenclature(CategoricalFusion(**fusion_example))
    assert nm == fixture_nomenclature

    nm = fusor.generate_nomenclature(CategoricalFusion(**exhaustive_example))
    assert nm == fixture_nomenclature

    with open(EXAMPLES_DIR / "alk.json", "r") as alk_file:
        fusion = CategoricalFusion(**json.load(alk_file))
    nm = fusor.generate_nomenclature(fusion)
    assert nm == "ALK(hgnc:427)::v"

    with open(EXAMPLES_DIR / "epcam_msh2.json", "r") as epcam_file:
        fusion = CategoricalFusion(**json.load(epcam_file))
    nm = fusor.generate_nomenclature(fusion)
    assert nm == "refseq:NM_002354.2(EPCAM):e._5::AGGCTCCCTTGG::refseq:NM_000251.2(MSH2):e.2_"  # noqa: E501

    with open(EXAMPLES_DIR / "tpm3_ntrk1.json", "r") as ntrk_file:
        fusion = AssayedFusion(**json.load(ntrk_file))
    nm = fusor.generate_nomenclature(fusion)
    assert nm == "refseq:NM_152263.3(TPM3):e._8::refseq:NM_002529.3(NTRK1):e.10_"  # noqa: E501

    with open(EXAMPLES_DIR / "tpm3_pdgfrb.json", "r") as pdgfrb_file:
        fusion = CategoricalFusion(**json.load(pdgfrb_file))
    nm = fusor.generate_nomenclature(fusion)
    assert nm == "refseq:NM_152263.3(TPM3):e._8::refseq:NM_002609.3(PDGFRB):e.11_"  # noqa: E501

    with open(EXAMPLES_DIR / "ewsr1.json", "r") as ewsr1_file:
        fusion = AssayedFusion(**json.load(ewsr1_file))
    nm = fusor.generate_nomenclature(fusion)
    assert nm == "EWSR1(hgnc:3508)::?"

    nm = fusor.generate_nomenclature(igh_myc)
    assert nm == "reg_promoter@IGH(hgnc:5477)::MYC(hgnc:7553)"
