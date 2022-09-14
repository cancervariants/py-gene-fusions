"""Test nomenclature generation."""
from fusor.models import CategoricalFusion


def test_generate_nomenclature(fusor_instance, fusion_example, exhaustive_example):
    """Test that nomenclature generation is correct."""
    fixture_nomenclature = "reg_p@BRAF(hgnc:1097)::refseq:NM_152263.3(TPM3):e.1_8::ALK(hgnc:427)::ACGT::refseq:NC_000023.11(chr X):g.44908820_44908822(+)::v"  # noqa: E501
    nm = fusor_instance.generate_nomenclature(CategoricalFusion(**fusion_example))
    assert nm == fixture_nomenclature
    nm = fusor_instance.generate_nomenclature(CategoricalFusion(**exhaustive_example))
    assert nm == fixture_nomenclature

    from fusor import examples

    nm = fusor_instance.generate_nomenclature(examples.bcr_abl1)
    assert nm == "refseq:NM_004327.3(BCR):e.2+182::ACTAAAGCG::refseq:NM_005157.5(ABL1):e.2-173"  # noqa: E501

    nm = fusor_instance.generate_nomenclature(examples.bcr_abl1_expanded)
    assert nm == "refseq:NM_004327.3(BCR):e.2+182::ACTAAAGCG::refseq:NM_005157.5(ABL1):e.2-173"  # noqa: E501

    nm = fusor_instance.generate_nomenclature(examples.alk)
    assert nm == "ALK(hgnc:427)::v"

    nm = fusor_instance.generate_nomenclature(examples.tpm3_ntrk1)
    assert nm == "refseq:NM_152263.3(TPM3):e.8::refseq:NM_002529.3(NTRK1):e.10"  # noqa: E501

    nm = fusor_instance.generate_nomenclature(examples.tpm3_pdgfrb)
    assert nm == "refseq:NM_152263.3(TPM3):e.8::refseq:NM_002609.3(PDGFRB):e.11"  # noqa: E501

    nm = fusor_instance.generate_nomenclature(examples.ewsr1)
    assert nm == "EWSR1(hgnc:3508)::?"

    nm = fusor_instance.generate_nomenclature(examples.igh_myc)
    assert nm == "reg_e_EH38E3121735@IGH(hgnc:5477)::MYC(hgnc:7553)"
