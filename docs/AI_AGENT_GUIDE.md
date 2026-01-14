# Coding Agent Guide — hoodini

This guide gives coding agents a fast, safe overview of hoodini: architecture, entry points, key modules, data flows, configuration, outputs, and high-impact maintainability actions.
## Overview

- Purpose: large-scale, gene-centric comparative genomics toolkit. Downloads public assemblies, extracts gene neighborhoods, runs protein and nucleotide comparisons, merges external annotations, and builds trees to interpret genomic context.
- Packaging: Python with `setuptools`; console script `hoodini` declared in `pyproject.toml`.
- Inputs: single-column accession list or TSV sheet; optional TOML config.
- Outputs: unified tables (GFF, protein metadata, tree metadata), Newick tree, protein/nt link tables, and optional extra annotations.
## Entry Points

- Console entry: [pyproject.toml](../pyproject.toml#L1-L40) → `hoodini = hoodini.__main__:main`.
- Module shim: [hoodini/__main__.py](../hoodini/__main__.py#L1-L20) imports and runs `hoodini.cli.main`.
- CLI implementation: [hoodini/cli.py](../hoodini/cli.py)
  - `hoodini run`: main pipeline (inputs, windowing, threads, tree modes, protein/nt links, extra annotations).
  - `hoodini download`: support databases (`assembly_summary`, `metacerberus`, `type_dive`, `contig_lengths`, `databases`).
  - `hoodini utils`: utilities (`nuc2asmlen`, `prefetch_links`).
- Additional scripts: `deep_mmseqs`, `deep_jackhmmer` registered in [pyproject.toml](../pyproject.toml#L1-L40).

## Pipeline Flow (`hoodini run`)
1. Initialize inputs and output folder: [hoodini/initialize.py](../hoodini/initialize.py) → `initialize_inputs()`.
2. Parse IPG records: [hoodini/parse_ipg.py](../hoodini/parse_ipg.py) → `run_ipg()`.
3. Download/parse assemblies + extract neighborhoods: [hoodini/parse_assemblies.py](../hoodini/parse_assemblies.py) → `run_assembly_parser()`.
4. Protein all-vs-all (optional): [hoodini/protein_links.py](../hoodini/protein_links.py) → `run_protein_links()`.
5. Nucleotide pairwise / ANI (optional): [hoodini/pairwise_nt.py](../hoodini/pairwise_nt.py) → `run_pairwise_nt()`.
6. Cluster neighbor proteins: [hoodini/cluster_proteins.py](../hoodini/cluster_proteins.py) → `cluster_proteins()`.
7. Proteome similarity (AAI): [hoodini/pipeline/proteome_similarity.py](../src/hoodini/pipeline/proteome_similarity.py) → `run_proteome_similarity()`.
8. Taxonomy + tree building: [hoodini/taxonomy.py](../hoodini/taxonomy.py) → `parse_taxonomy_and_build_tree()`.
9. Extra annotations (optional): wrappers in [hoodini/extra_tools](../hoodini/extra_tools) for PADLOC, DefenseFinder, eggNOG-mapper (emapper), GenoMAD, ncRNA/Infernal, CCtyper, domains, PHROGs, anti-defense.
10. Write unified outputs for visualization: [hoodini/write_data.py](../hoodini/write_data.py) → `write_viz_outputs()`.

Execution notes:
- Effective config precedence: `config/defaults.toml` < user TOML (`--config`) < explicit CLI flags.
- Flags gate stages, e.g., `--tree-mode aai_tree` enables AAI; `--nt-links`/`--ani-mode` enable NT/ANI.

## Configuration
- Defaults: [hoodini/config/defaults.toml](../hoodini/config/defaults.toml) loaded via `hoodini.config.load_default_config()`.
- Schema/validation: [hoodini/config/schema.py](../hoodini/config/schema.py).
- Examples: [hoodini/config/example.toml](../hoodini/config/example.toml), [hoodini/config/metacerberus.toml](../hoodini/config/metacerberus.toml).
- Merge logic in CLI: see [hoodini/cli.py](../hoodini/cli.py#L1-L200), uses `tomli` and `click.ParameterSource` to keep only user-set flags.

## Key Modules
- CLI + helpers: [hoodini/cli.py](../hoodini/cli.py), [hoodini/utils/cli_helpers.py](../hoodini/utils/cli_helpers.py).
- Core utilities: [hoodini/utils/core.py](../hoodini/utils/core.py), [hoodini/utils/file_formats.py](../hoodini/utils/file_formats.py), [hoodini/utils/logging_utils.py](../hoodini/utils/logging_utils.py), [hoodini/utils/ncbi_api.py](../hoodini/utils/ncbi_api.py).
- Downloads: [hoodini/download](../hoodini/download) — `assembly_summary.py`, `contig_lengths.py`, `metacerberus.py`, `type_dive.py`, `databases.py`.
- Packaged data: [hoodini/data](../hoodini/data) — parquets, emapper/mmseqs references, `genomad_db`, `all.cm`.
- External tool wrappers: [hoodini/extra_tools](../hoodini/extra_tools) — PADLOC, DefenseFinder, emapper, GenoMAD, BLAST, domains, CCtyper, PHROGs, ncRNA.
- Basic/extra plotting: [hoodini/basic_plot.py](../hoodini/basic_plot.py), [hoodini/extra_plot.py](../hoodini/extra_plot.py); template in [hoodini/template/template.html](../hoodini/template/template.html).

## CLI Utilities (`hoodini utils`)

- `nuc2asmlen`: [hoodini/pipeline/helpers/nuc2asmlen.py](../src/hoodini/pipeline/helpers/nuc2asmlen.py) → assembly/contig length metadata for nuccore/contig accessions.
- `prefetch_links`: [hoodini/pipeline/helpers/prefetch_links.py](../src/hoodini/pipeline/helpers/prefetch_links.py) → table of `assembly_id`, `file_type`, `link`.

## Downloads (`hoodini download`)

- `assembly_summary`: [hoodini/download/assembly_summary.py](../hoodini/download/assembly_summary.py) → local parquet.
- `contig_lengths`: [hoodini/download/contig_lengths.py](../hoodini/download/contig_lengths.py).
- `metacerberus`: [hoodini/download/metacerberus.py](../hoodini/download/metacerberus.py).
- `type_dive`: [hoodini/download/type_dive.py](../hoodini/download/type_dive.py).
- `databases`: [hoodini/download/databases.py](../hoodini/download/databases.py) → emapper/mmseqs DB, PADLOC models, DefenseFinder models, GenoMAD DB, eggNOG parquet support.

## Standard Outputs
Visualization files are written under `results/hoodini-viz/` with this structure:

- Root:
  - `tree.nwk` — Newick-formatted tree string.
  - `hoodini-viz.html` — standalone viewer (embeds base64-compressed parquet data).
- `tsv/` (tabular text):
  - `gff.gff` — combined GFF for parsed assemblies.
  - `hoods.txt` — hoods per neighborhood (`hood_id`, `seqid`, `start`, `end`, `align_gene`).
  - `protein_metadata.txt` — protein metadata (clusters, product, merged annotations).
  - `tree_metadata.txt` — per-leaf tree metadata.
  - `nucleotide_links.txt` — nucleotide links; placeholder if not produced.
  - `protein_links.txt` — protein links; placeholder if not produced.
- `parquet/` (columnar):
  - `gff.parquet`, `hoods.parquet`, `protein_metadata.parquet`, `tree_metadata.parquet`, `nucleotide_links.parquet`, `protein_links.parquet`, and optional `domains.parquet`, `domains_metadata.parquet` when domains are requested.

## Dependencies & Environment
- Conda/Mamba: see [environment.yml](../environment.yml). Includes `bioconda` tools (PADLOC, GenoMAD, MAFFT, Foldseek, FastANI, VeryFastTree, Infernal, FAMSA, DIAMOND, Skani) and pip packages (pyhmmer, diamondonpy, aria2p, jinja2, etc.).
- Python packages: pandas, numpy, matplotlib, rich, rich-click, pyarrow, polars, biopython, networkx, requests, tomli, ete3, httpx, aiohttp.

## Development Standards (Recommended)
- Typing: add type hints and `mypy` to public functions, especially under `utils/` and `extra_tools/`.
- Lint/format: integrate `ruff` + `black` + `isort` and set up pre-commit hooks.
- Testing: add `pytest` with small fixtures; unit tests for parsers and runners (`initialize`, `parse_ipg`, `parse_assemblies`, `pairwise_nt`, `cluster_proteins`, `taxonomy`, `write_data`).
- Logging: standardize levels/messages via `logging_utils`; avoid direct prints.
- Errors: raise specific exceptions with actionable messages; validate inputs early.
- Config: consolidate validation in `schema.py`; document keys in `defaults.toml` and `README`.
- IO contracts: prefer `pyarrow`/`polars` for parquet; enforce consistent schemas.

## Refactor Opportunities (Priority)
1. CLI: move heavy logic out of [hoodini/cli.py](../hoodini/cli.py) into dedicated pipeline/services modules to reduce coupling and ease testing.
2. Config model: replace ad-hoc `ConfigObj` with `@dataclass` (or pydantic), with explicit validation and separated merge/resolution.
3. `extra_tools`: normalize `run_*` interfaces (consistent signature: inputs, `output_dir`, `threads`, `uids`) and error/return handling.
4. `pairwise_nt` / `protein_links`: decouple engines (blastn/fastANI/minimap2) via strategy pattern; add deterministic tests with mini datasets.
5. `write_data`: centralize column schemas and validators; document input/output contract.
6. `download/*`: factor shared code (progress, retries, checksum) and parameterize destinations.
7. Documentation: auto-generate API docs (pdoc/sphinx) and reproducible examples under [example/](../example).

## Quick Start for Agents
- Read first: [README.md](../README.md), [pyproject.toml](../pyproject.toml), [environment.yml](../environment.yml), [hoodini/cli.py](../hoodini/cli.py).
- Typical modifications:
  - Add a new annotation: create `hoodini/extra_tools/new_tool.py` with `run_new_tool(...)` and wire into `cli.run`.
  - Change parameter resolution: edit `hoodini/config/schema.py` and `defaults.toml`; update merge in `cli.py`.
  - Extend comparison engines: add a strategy to `pairwise_nt.py` or `protein_links.py`.
- Quick tests:
  - Use [example/input.txt](../example/input.txt) and [example/hoodini_cli_workflow.ipynb](../example/hoodini_cli_workflow.ipynb).
  - Minimal CLI: `hoodini run --input example/input.txt --output example/results --num-threads 2 --force`.

## Risks & Assumptions
- External tools (PADLOC, GenoMAD, emapper, etc.) must be installed and available in the environment; some runners expect/download local DBs under `hoodini/data`.
- Network-heavy steps (NCBI/IMG) may be rate-limited; use `--api-key` and tune concurrency.
- Parquet/table schemas must remain stable; changes require migration and documentation updates.

## Maintainability Checklist (Actionable)
- Add `ruff`, `black`, `isort`, `mypy`, `pytest`, and `pre-commit`.
- Create `CONTRIBUTING.md` with style, commit conventions, and testing guide.
- Add CI (GitHub Actions) for lint, type-check, unit tests, and package build.
- Type public functions and ensure docstrings document signatures and side effects.
- Stabilize output schemas with clear contracts in `write_data`.
- Unify configuration as a `@dataclass` with centralized validation.

---

This document is designed for quick reading by coding agents: it links key locations and contracts without repeating source code. For deeper detail, follow the linked modules and examples in `example/`.
