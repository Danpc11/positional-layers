"""Pipeline behaviour that the analysis scripts rely on: idempotent result tables and an end-to-end run of the CRISPRi
model on synthetic pairs covering every distance class."""
import os, subprocess, sys, pathlib
import numpy as np, pandas as pd
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]


def test_upsert_replaces_rows_instead_of_appending(tmp_path, monkeypatch):
    monkeypatch.setenv('POSLAYERS_RESULTS', str(tmp_path))
    from poslayers.config import upsert_csv
    p = str(tmp_path / 't.csv')
    upsert_csv(pd.DataFrame({'tissue': ['a', 'b'], 'v': [1, 2]}), p, ['tissue'])
    upsert_csv(pd.DataFrame({'tissue': ['b'], 'v': [20]}), p, ['tissue'])          # re-run of tissue b
    out = pd.read_csv(p).set_index('tissue').v.to_dict()
    assert out == {'a': 1, 'b': 20}


def test_upsert_rejects_duplicate_keys(tmp_path):
    from poslayers.config import upsert_csv
    with pytest.raises(ValueError):
        upsert_csv(pd.DataFrame({'k': [1, 1], 'v': [1, 2]}), str(tmp_path / 'd.csv'), ['k'])


def test_crispri_analysis_runs_end_to_end(tmp_path):
    rng = np.random.default_rng(0); rows = []
    for p in range(60):
        kd = rng.uniform(1, 4)
        for typ in ('cis', 'trans'):
            n = 40
            d = rng.choice([5e3, 3e4, 1e5, 3e5, 8e5], n) if typ == 'cis' else np.full(n, np.nan)
            r = rng.uniform(-0.5, 0.8, n)
            resp = -1.5 * (np.nan_to_num(d, nan=1e9) < 1e4) * kd / 3 + rng.normal(0, 1, n)
            rows.append(pd.DataFrame({'pert': p, 'target': f'G{p}', 'type': typ, 'kd': kd, 'r': r, 'resp': resp, 'dist': d}))
    pd.concat(rows).to_csv(tmp_path / 'crispri_pairs.csv.gz', index=False)
    env = dict(os.environ, POSLAYERS_RESULTS=str(tmp_path), PYTHONPATH=str(ROOT / 'src'))
    res = subprocess.run([sys.executable, str(ROOT / 'scripts' / 'review' / 'crispri_analyse.py')], env=env, capture_output=True, text=True)
    assert res.returncode == 0, res.stderr[-2000:]
    assert 'coupling x knockdown' in res.stdout
    assert (tmp_path / 'crispri_coupling_distribution_cis_vs_trans.csv').exists()
