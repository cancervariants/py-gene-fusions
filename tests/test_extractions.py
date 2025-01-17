"""Module for testing extraction methods"""

from pathlib import Path

from fusor.extract import (
    get_arriba_records,
    get_cicero_records,
    get_fusion_catcher_records,
    get_jaffa_records,
    get_star_fusion_records,
)


def test_get_jaffa_records(fixture_data_dir):
    """Test that get_jaffa_records works correctly"""
    path = fixture_data_dir / "jaffa_results.csv"
    fusions_list = get_jaffa_records(Path(path))
    assert len(fusions_list) == 491

    path = fixture_data_dir / "jaffa_resultss.csv"
    fusions_list = get_jaffa_records(Path(path))
    assert fusions_list is None


def test_get_star_fusion_records(fixture_data_dir):
    """Test that get_star_fusion_records works correctly"""
    path = fixture_data_dir / "star-fusion.fusion_predictions.abridged.tsv"
    fusions_list = get_star_fusion_records(Path(path))
    assert len(fusions_list) == 37

    path = fixture_data_dir / "star-fusion.fusion_predictions.abridged.tsvs"
    fusions_list = get_star_fusion_records(Path(path))
    assert fusions_list is None


def test_get_fusion_catcher_records(fixture_data_dir):
    """Test that get_fusion_catcher_records works correctly"""
    path = fixture_data_dir / "final-list_candidate-fusion-genes.txt"
    fusions_list = get_fusion_catcher_records(Path(path))
    assert len(fusions_list) == 355

    path = fixture_data_dir / "final-list_candidate-fusion-genes.txts"
    fusions_list = get_fusion_catcher_records(Path(path))
    assert fusions_list is None


def test_get_arriba_records(fixture_data_dir):
    """Test that get_arriba_records works correctly"""
    path = fixture_data_dir / "fusions_arriba_test.tsv"
    fusions_list = get_arriba_records(Path(path))
    assert len(fusions_list) == 1

    path = fixture_data_dir / "fusionsd_arriba_test.tsv"
    fusions_list = get_arriba_records(Path(path))
    assert fusions_list is None


def test_get_cicero_records(fixture_data_dir):
    """Test that get_cicero_records works correctly"""
    path = fixture_data_dir / "annotated.fusion.txt"
    fusions_list = get_cicero_records(Path(path))
    assert len(fusions_list) == 1

    path = fixture_data_dir / "annnotated.fusion.txt"
    fusions_list = get_cicero_records(Path(path))
    assert fusions_list is None
