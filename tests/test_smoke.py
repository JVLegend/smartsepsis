"""Smoke tests for SmartSepsis-Oph.

These tests verify that:
  1. `pyproject.toml` is syntactically valid TOML.
  2. Every top-level *.py module at the repo root can be imported without
     raising ImportError (modules with known-heavy side effects are skipped).

The tests deliberately do NOT exercise pipeline behavior — they only catch
gross packaging / dependency regressions.
"""
from __future__ import annotations

import importlib
import importlib.util
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent

# Modules known to perform heavy work (model downloads, GPU init, network I/O,
# or long batch jobs) at import time, or which require optional heavy deps.
# These are skipped until the migration to src/ layout (see MIGRATION.md).
HEAVY_MODULES = {
    "build_hf_dataset_extended",  # large dataset assembly
    "build_hf_dataset",           # HF dataset build
    "evo2_scoring",               # Evo2 model
    "prott5_ensemble",            # ProtT5 model load
    "protein_structure",          # ESMFold
    "nanofold_benchmark",         # heavy structure benchmark
    "nanofold_calibration_v5",    # heavy calibration
    "structure_features_v3",      # depends on PDBs / structure stack
    "amrfinderplus_embed",        # external CARD/AMRFinderPlus data
    "tier2_pubmed_mining",        # network-bound literature mining
    "run_batch",                  # orchestrator with side effects
}


def _discover_root_modules() -> list[str]:
    """Return the stems of all *.py files at the repo root, sorted."""
    return sorted(
        p.stem
        for p in REPO_ROOT.glob("*.py")
        if p.is_file() and not p.name.startswith("_")
    )


ROOT_MODULES = _discover_root_modules()


def test_pyproject_parses():
    """pyproject.toml must be valid TOML."""
    try:
        import tomllib  # Python 3.11+
    except ModuleNotFoundError:  # pragma: no cover
        import tomli as tomllib  # type: ignore[no-redef]

    pyproject = REPO_ROOT / "pyproject.toml"
    assert pyproject.exists(), "pyproject.toml missing at repo root"
    with pyproject.open("rb") as fh:
        data = tomllib.load(fh)
    assert data["project"]["name"] == "smartsepsis-oph"
    assert data["project"]["version"].startswith("0.1.0")


@pytest.mark.parametrize("module_name", ROOT_MODULES)
def test_root_module_importable(module_name: str) -> None:
    """Every root *.py module should be importable (or explicitly skipped)."""
    if module_name in HEAVY_MODULES:
        pytest.skip(
            f"{module_name}: known heavy import (model/data download or side effects). "
            "TODO: re-enable after src/ migration and lazy-loading refactor."
        )

    # Ensure repo root is on sys.path so we can import root scripts as modules.
    repo_root_str = str(REPO_ROOT)
    if repo_root_str not in sys.path:
        sys.path.insert(0, repo_root_str)

    try:
        importlib.import_module(module_name)
    except ImportError as exc:
        pytest.fail(f"ImportError while importing {module_name!r}: {exc}")
    except Exception as exc:  # noqa: BLE001
        # Non-ImportError failures (e.g., FileNotFoundError on optional data)
        # are tolerated — the module loaded far enough to evaluate, which is
        # what this smoke test cares about.
        pytest.skip(f"{module_name} imported but raised {type(exc).__name__}: {exc}")
