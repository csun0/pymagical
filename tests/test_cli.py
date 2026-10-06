"""CLI diagnostics for input paths, without running inference."""
import sys
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from pymagical.cli import main


@pytest.fixture
def inference_mock(monkeypatch):
    run = Mock()
    monkeypatch.setitem(sys.modules, 'pymagical.magical', SimpleNamespace(run_magical=run))
    return run


def test_independent_cell_directory_reports_only_missing_paths(
    tmp_path, monkeypatch, capsys, inference_mock
):
    monkeypatch.chdir(tmp_path)
    data = tmp_path / 'Demo_input_files'
    data.mkdir()
    for filename in (
        'Motif_mapping_prior.txt', 'Motifs.txt',
        'RaoGM12878_40kb_TopDomTADs_filtered_hg38.txt', 'hg38_Refseq.txt',
    ):
        (data / filename).touch()
    monkeypatch.setattr(sys, 'argv', [
        'pymagical', 'run', '--main-dir', 'Demo_input_files/',
        '--cell-dir', 'Demo_input_files/',
    ])
    with pytest.raises(SystemExit) as error:
        main()
    assert error.value.code == 1
    output = capsys.readouterr().out
    assert '8 required input file(s) missing or invalid' in output
    expected = data / 'Cell_type_candidate_genes.txt'
    assert f'[MISSING] Candidate Genes (--cand-genes): {expected}' in output
    assert 'Motif Info' not in output
    assert 'Motif Mapping' not in output
    assert 'RefSeq File' not in output
    assert 'TAD File' not in output
    assert '[FOUND]' not in output
    assert f'Resolved cell directory: {data}' in output
    assert '--main-dir and --cell-dir are independent paths' in output
    inference_mock.assert_not_called()


def test_unset_inputs_show_override_flags(monkeypatch, capsys, inference_mock):
    monkeypatch.setattr(sys, 'argv', ['pymagical', 'run', '--no-tad'])
    with pytest.raises(SystemExit):
        main()
    output = capsys.readouterr().out
    assert '[UNSET] RNA Counts (--rna-counts): no path supplied' in output
    assert 'TAD File' not in output
    inference_mock.assert_not_called()


def test_directory_override_is_rejected(tmp_path, monkeypatch, capsys, inference_mock):
    monkeypatch.setattr(sys, 'argv', [
        'pymagical', 'run', '--cand-genes', str(tmp_path),
    ])
    with pytest.raises(SystemExit):
        main()
    assert f'[NOT A FILE] Candidate Genes (--cand-genes): {tmp_path}' in capsys.readouterr().out
    inference_mock.assert_not_called()


@pytest.mark.parametrize('absolute', [False, True])
def test_separate_directory_paths_reach_inference(
    tmp_path, monkeypatch, inference_mock, absolute
):
    monkeypatch.chdir(tmp_path)
    shared = tmp_path / 'shared'
    cells = tmp_path / 'cells'
    shared.mkdir()
    cells.mkdir()
    for filename in (
        'Motif_mapping_prior.txt', 'Motifs.txt',
        'RaoGM12878_40kb_TopDomTADs_filtered_hg38.txt', 'hg38_Refseq.txt',
    ):
        (shared / filename).touch()
    for filename in (
        'Cell_type_candidate_genes.txt', 'Cell_type_candidate_peaks.txt',
        'Cell_type_scRNA_read_count.txt', 'scRNA_genes.txt',
        'Cell_type_scRNA_cell_meta.txt', 'Cell_type_scATAC_read_count.txt',
        'scATAC_peaks.txt', 'Cell_type_scATAC_cell_meta.txt',
    ):
        (cells / filename).touch()
    monkeypatch.setattr(sys, 'argv', [
        'pymagical', 'run', '--main-dir', str(shared) if absolute else 'shared',
        '--cell-dir', str(cells) if absolute else 'cells',
    ])
    main()
    inference_mock.assert_called_once()
    from pathlib import Path
    arguments = inference_mock.call_args.kwargs
    assert Path(arguments['rna_counts_file']).resolve() == cells / 'Cell_type_scRNA_read_count.txt'
    assert Path(arguments['motif_name_file']).resolve() == shared / 'Motifs.txt'
