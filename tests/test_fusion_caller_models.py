"""Module for testing fusion caller classes"""

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
    records = JAFFA.load_records(path)
    assert len(records) == 491

    path = Path(fixture_data_dir / "jaffa_resultss.csv")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert JAFFA.load_records(path)


def test_get_star_fusion_records(fixture_data_dir):
    """Test that get_star_fusion_records works correctly"""
    path = Path(fixture_data_dir / "star-fusion.fusion_predictions.abridged.tsv")
    records = STARFusion.load_records(path)
    assert len(records) == 37

    path = Path(fixture_data_dir / "star-fusion.fusion_predictions.abridged.tsvs")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert STARFusion.load_records(path)


def test_get_fusion_catcher_records(fixture_data_dir):
    """Test that get_fusion_catcher_records works correctly"""
    path = Path(fixture_data_dir / "final-list_candidate-fusion-genes.txt")
    fusions_list = FusionCatcher.load_records(path)
    assert len(fusions_list) == 355

    path = Path(fixture_data_dir / "final-list_candidate-fusion-genes.txts")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert FusionCatcher.load_records(path)


def test_get_arriba_records(fixture_data_dir):
    """Test that get_arriba_records works correctly"""
    path = Path(fixture_data_dir / "fusions_arriba_test.tsv")
    fusions_list = Arriba.load_records(path)
    assert len(fusions_list) == 1

    path = Path(fixture_data_dir / "fusionsd_arriba_test.tsv")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert Arriba.load_records(path)


def test_get_cicero_records(fixture_data_dir):
    """Test that get_cicero_records works correctly"""
    path = Path(fixture_data_dir / "annotated.fusion.txt")
    fusions_list = Cicero.load_records(path)
    assert len(fusions_list) == 1

    path = Path(fixture_data_dir / "annnotated.fusion.txt")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert Cicero.load_records(path)


def test_get_enfusion_records(fixture_data_dir):
    """Test that get_enfusion_records works correctly"""
    path = Path(fixture_data_dir / "enfusion_test.csv")
    fusions_list = EnFusion.load_records(path)
    assert len(fusions_list) == 1

    path = Path(fixture_data_dir / "enfusions_test.csv")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert EnFusion.load_records(path)


def test_get_genie_records(fixture_data_dir):
    """Test that get_genie_records works correctly"""
    path = Path(fixture_data_dir / "genie_test.txt")
    fusions_list = Genie.load_records(path)
    assert len(fusions_list) == 1

    path = Path(fixture_data_dir / "genie_tests.txt")
    with pytest.raises(ValueError, match=f"{path} does not exist"):
        assert Genie.load_records(path)
