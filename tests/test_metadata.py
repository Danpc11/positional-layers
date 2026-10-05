"""Release metadata must agree: package, pyproject, CITATION.cff and .zenodo.json carry the same version."""
import json, pathlib, re
import poslayers

ROOT = pathlib.Path(__file__).resolve().parents[1]


def test_versions_agree():
    v = poslayers.__version__
    assert re.search(r'^version = "([^"]+)"', (ROOT / "pyproject.toml").read_text(), re.M).group(1) == v
    assert re.search(r"^version: ['\"]?([^'\"\n]+)", (ROOT / "CITATION.cff").read_text(), re.M).group(1).strip() == v
    assert json.loads((ROOT / ".zenodo.json").read_text())["version"] == v


def test_readme_version_badge_matches():
    m = re.search(r"img\.shields\.io/badge/version-([0-9.]+)-", (ROOT / "README.md").read_text())
    assert m and m.group(1) == poslayers.__version__
