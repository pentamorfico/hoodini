# Hoodini Architecture and Refactor Plan

This document captures the target structure and conventions to standardize the project. The
package now lives under `src/` (PEP 517/518 layout); future moves should keep that root.

## Target Directory Layout
- `src/hoodini/` (move package under `src/` for clean imports)
  - `cli/`: click entrypoints only (parse args, build config, call pipeline)
  - `config/`: defaults, schema, validation, settings
  - `pipeline/`: orchestrated stages (`inputs`, `ipg`, `assemblies`, `comparisons`, `clustering`, `taxonomy`, `annotations`, `outputs`)
  - `services/`: external integrations (NCBI/IMG fetchers, tool runners, downloaders)
  - `data_access/`: readers/writers for parquet/fasta/gff
  - `models/`: dataclasses/enums/Polars schemas
  - `utils/`: focused helpers (io, validation, logging, polars adapters)
  - `templates/`, `assets/`: bundled HTML/figures if needed
- `docs/`, `tests/`, `example/`, `scripts/` remain at repo root.
- Bundled data: prefer download-on-demand; keep only minimal fixtures in `tests/data`.

## Configuration Flow
1) Load defaults from `config/defaults.toml`.
2) Merge user TOML (optional) and explicit CLI flags into `RuntimeConfig`.
3) Pass typed config into pipeline stages; avoid ad-hoc dicts.

## Pipeline Responsibilities
- `inputs`: validate/normalize IDs, create working directory (no prompts inside library code).
- `ipg`: fetch/parse IPG and resolve representatives.
- `assemblies`: download/parse assemblies, extract neighborhoods.
- `comparisons`: protein and nucleotide pairwise computations.
- `clustering`: protein family clustering, sORF handling.
- `taxonomy`: tree building and metadata enrichment.
- `annotations`: domain/blast/padloc/defensefinder/emapper/genomad/etc.
- `outputs`: write combined viz-ready tables and trees.

Each stage should accept/return typed dataclasses (or Polars DataFrames validated against `models/schemas.py`).

## Naming & API Surface
- Public API modules should avoid leading underscores; keep private helpers `_prefixed`.
- Reserve `run_*` for CLI-facing entry wrappers; internal helpers should use descriptive verbs.
- Export a stable surface in `hoodini/__init__.py` for integrations; keep CLI UX in `hoodini/cli`.

## Testing Strategy
- Unit tests per pipeline stage (small Polars fixtures).
- Integration smoke test guarded by env flag (already present).
- Add fast tests for config merge logic and CLI validation.

## Near-Term Refactor Checklist
- Move package to `src/` layout and update `pyproject.toml` accordingly.
- Split `hoodini/cli.py` into `cli/main.py` (group/command definitions) and `pipeline/runner.py` (orchestration).
- Break up `hoodini/utils/core.py` into smaller modules: `validation.py`, `id_parsing.py`, `io.py`, `logging_setup.py`.
- Remove import-time side effects and logging configuration from libraries; configure logging in CLI only.
- Align dependencies with Python 3.10+ (matches black/ruff config) and trim unused packages from `environment.yml`.
