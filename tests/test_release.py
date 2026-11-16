"""DepMap release provenance checks (README header, RELEASE.txt, Model.csv schema)."""
from types import SimpleNamespace

from src.phase0 import download, load

README = """# DepMap Public 26Q3

### Model.csv

- ModelID: id
- PatientID: patient

### OmicsExpression.csv
- Foo: bar
"""


def _setup(tmp_path, monkeypatch, readme=README, release="DepMap Public 26Q3, downloaded 2026-11-16",
           header="ModelID,PatientID\n"):
    (tmp_path / "README.txt").write_text(readme)
    (tmp_path / "Model.csv").write_text(header)
    rel = tmp_path / "RELEASE.txt"
    rel.write_text(release)
    srcs = {"model": SimpleNamespace(path=tmp_path / "Model.csv"),
            "model_readme": SimpleNamespace(path=tmp_path / "README.txt")}
    monkeypatch.setattr(load, "SOURCES", srcs)
    monkeypatch.setattr(load, "MODEL_RELEASE_FILE", rel)
    monkeypatch.setattr(download, "MODEL_RELEASE_FILE", rel)
    monkeypatch.setattr(download, "REPO", tmp_path)


def test_consistent(tmp_path, monkeypatch):
    _setup(tmp_path, monkeypatch)
    ev = load.release_evidence()
    assert ev["readme_release"] == "DepMap Public 26Q3" and ev["schema_match"] is True
    assert download.release_problems() == []
    assert load.model_release().startswith("DepMap Public 26Q3 (from depmap/README.txt header")


def test_release_txt_disagrees(tmp_path, monkeypatch):
    _setup(tmp_path, monkeypatch, release="DepMap Public 26Q1")
    assert any("disagrees" in p for p in download.release_problems())


def test_schema_mismatch(tmp_path, monkeypatch):
    _setup(tmp_path, monkeypatch, header="ModelID,PatientID,NewColumn\n")
    assert any("do not match" in p for p in download.release_problems())


def test_drop_par_y():
    import pandas as pd
    df = pd.DataFrame({"a": [1.0, 0.0, 2.0]}, index=["ENSG00000000003.10", "ENSGR0000002586.13", "ENSG00000002586.13"])
    assert list(load.drop_par_y(df).index) == ["ENSG00000000003.10", "ENSG00000002586.13"]
