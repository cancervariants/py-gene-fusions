"""Module for testing extraction methods"""

from pathlib import Path

import pytest

from fusor.fusion_caller_models import (
    JAFFA,
    Arriba,
    Cicero,
    EnFusion,
    FusionCatcher,
    Genie,
    STARFusion,
)


def test_get_jaffa_records(fixture_data_dir):
    """Test that get_jaffa_records works correctly"""
    path = Path(fixture_data_dir / "jaffa_results.csv")
    jaffa_instance = JAFFA()
    records = jaffa_instance.load_records(path)
    assert len(records) == 491

    path = Path(fixture_data_dir / "jaffa_resultss.csv")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert jaffa_instance.load_records(path)


def test_get_star_fusion_records(fixture_data_dir):
    """Test that get_star_fusion_records works correctly"""
    path = Path(fixture_data_dir / "star-fusion.fusion_predictions.abridged.tsv")
    sf_instance = STARFusion()
    records = sf_instance.load_records(path)
    assert len(records) == 37

    path = Path(fixture_data_dir / "star-fusion.fusion_predictions.abridged.tsvs")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert sf_instance.load_records(path)

    # path = fixture_data_dir / "star-fusion.fusion_predictions.abridged.tsv"
    # fusions_list = get_star_fusion_records(Path(path))
    # assert len(fusions_list) == 37

    # path = fixture_data_dir / "star-fusion.fusion_predictions.abridged.tsvs"
    # with pytest.raises(ValueError, match=f"{path} does not exist"):
    #   assert get_jaffa_records(path)


def test_get_fusion_catcher_records(fixture_data_dir):
    """Test that get_fusion_catcher_records works correctly"""
    path = Path(fixture_data_dir / "final-list_candidate-fusion-genes.txt")
    fc_instance = FusionCatcher()
    fusions_list = fc_instance.load_records(path)
    assert len(fusions_list) == 355

    path = Path(fixture_data_dir / "final-list_candidate-fusion-genes.txts")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert fc_instance.load_records(path)


def test_get_arriba_records(fixture_data_dir):
    """Test that get_arriba_records works correctly"""
    path = Path(fixture_data_dir / "fusions_arriba_test.tsv")
    arriba = Arriba()
    fusions_list = arriba.load_records(path)
    assert len(fusions_list) == 1

    path = Path(fixture_data_dir / "fusionsd_arriba_test.tsv")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert arriba.load_records(path)


def test_get_cicero_records(fixture_data_dir):
    """Test that get_cicero_records works correctly"""
    path = Path(fixture_data_dir / "annotated.fusion.txt")
    cicero = Cicero()
    fusions_list = cicero.load_records(path)
    assert len(fusions_list) == 1

    path = Path(fixture_data_dir / "annnotated.fusion.txt")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert cicero.load_records(path)


def test_get_enfusion_records(fixture_data_dir):
    """Test that get_enfusion_records works correctly"""
    path = Path(fixture_data_dir / "enfusion_test.csv")
    enfusion = EnFusion()
    fusions_list = enfusion.load_records(path)
    assert len(fusions_list) == 1

    path = Path(fixture_data_dir / "enfusions_test.csv")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert enfusion.load_records(path)


def test_get_genie_records(fixture_data_dir):
    """Test that get_genie_records works correctly"""
    path = Path(fixture_data_dir / "genie_test.txt")
    genie = Genie()
    fusions_list = genie.load_records(path)
    assert len(fusions_list) == 1

    path = Path(fixture_data_dir / "genie_tests.txt")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert genie.load_records(path)
