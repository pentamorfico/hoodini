# AGENTS.md — Hoodini

> Context file for AI coding agents (Copilot, Cursor, Claude, etc.)

## Setup Commands

```bash
# Create and activate environment (REQUIRED)
mamba env create -f environment.yml
mamba activate hoodini

# Install in editable mode for development
pip install -e ".[dev]"

# Run tests
pytest tests/

# Lint and format
black . && isort . && ruff check .
```

> ⚠️ **Important:** Always activate the conda environment before running any commands. External bioinformatics tools (PADLOC, MMseqs2, etc.) won't work without it.

## Project Overview

**Hoodini** is a large-scale, gene-centric comparative genomics toolkit. It downloads public assemblies, extracts gene neighborhoods, runs protein and nucleotide comparisons, merges external annotations, and builds trees to interpret genomic context.

- **Language:** Python 3.10+
- **Packaging:** `setuptools` with console script `hoodini` declared in `pyproject.toml`
- **Inputs:** Single-column accession list or TSV sheet; optional TOML config
- **Outputs:** Unified tables (GFF, protein metadata, tree metadata), Newick tree, protein/nt link tables, and optional extra annotations

## Entry Points

| Entry | Location | Purpose |
|-------|----------|---------|
| Console | `pyproject.toml` → `hoodini = hoodini.__main__:main` | Main CLI entry |
| Module shim | `hoodini/__main__.py` | Imports and runs `hoodini.cli.main` |
| CLI | `hoodini/cli.py` | Command definitions |

### CLI Commands

```bash
hoodini run        # Main pipeline
hoodini download   # Download support databases
hoodini utils      # Utility commands (nuc2asmlen, prefetch_links)
```

## Architecture

```
hoodini/
├── cli.py                 # Click entrypoints (parse args, build config, call pipeline)
├── config/                # defaults.toml, schema.py, validation
├── pipeline/              # Orchestrated stages
├── download/              # Database downloaders
├── extra_tools/           # External tool wrappers (PADLOC, DefenseFinder, etc.)
├── utils/                 # Helpers (io, validation, logging, NCBI API)
├── data/                  # Bundled parquets, references
└── template/              # HTML template for visualization
```

## Pipeline Flow (`hoodini run`)

```
1. initialize.py      → Initialize inputs and output folder
2. parse_ipg.py       → Parse IPG records
3. parse_assemblies.py → Download/parse assemblies, extract neighborhoods
4. protein_links.py   → Protein all-vs-all (optional)
5. pairwise_nt.py     → Nucleotide pairwise / ANI (optional)
6. cluster_proteins.py → Cluster neighbor proteins
7. proteome_similarity.py → AAI computation
8. taxonomy.py        → Taxonomy + tree building
9. extra_tools/*      → Optional annotations (PADLOC, DefenseFinder, etc.)
10. write_data.py     → Write unified outputs for visualization
```

## Configuration

- **Defaults:** `hoodini/config/defaults.toml`
- **Schema:** `hoodini/config/schema.py`
- **Precedence:** `defaults.toml` < user TOML (`--config`) < CLI flags

## Key Modules

| Module | Purpose |
|--------|---------|
| `cli.py` | CLI definitions and orchestration |
| `utils/core.py` | Core utilities |
| `utils/ncbi_api.py` | NCBI API interactions |
| `download/*` | Database downloaders |
| `extra_tools/*` | External tool wrappers |

## Output Structure

```
results/hoodini-viz/
├── tree.nwk                    # Newick tree
├── hoodini-viz.html            # Standalone viewer
├── tsv/
│   ├── gff.gff
│   ├── hoods.txt
│   ├── protein_metadata.txt
│   ├── tree_metadata.txt
│   ├── nucleotide_links.txt
│   └── protein_links.txt
└── parquet/
    └── (same files as .parquet)
```

## Environment & Dependencies

> ⚠️ **Important:** Hoodini requires a conda/mamba environment to run correctly. Many features depend on external bioinformatics tools that must be installed via Bioconda.

```bash
# Create and activate environment
mamba env create -f environment.yml
mamba activate hoodini
```

- **Conda/Mamba:** See `environment.yml` for full dependency list
- **Bioconda tools:** PADLOC, GenoMAD, MAFFT, Foldseek, FastANI, VeryFastTree, Infernal, FAMSA, DIAMOND, Skani, MMseqs2, pyhmmer
- **Python packages:** pandas, polars, pyarrow, biopython, rich, click, httpx, ete3

Running without the environment will cause errors when pipeline stages try to call external tools.

## Testing

```bash
# Run all tests
pytest tests/

# Run with verbose output
pytest tests/ -v --tb=short

# Run specific test file
pytest tests/test_config.py

# Run integration tests (slower, requires databases)
HOODINI_RUN_SLOW_TESTS=1 pytest tests/
```

## Development Priorities

1. **CLI refactor:** Move heavy logic out of `cli.py` into dedicated pipeline modules
2. **Config model:** Replace ad-hoc dicts with typed dataclasses
3. **Testing:** Add pytest fixtures for pipeline stages
4. **extra_tools:** Normalize `run_*` interfaces across all wrappers
5. **Documentation:** Auto-generate API docs

## Common Tasks

### Adding a new annotation tool

1. Create wrapper in `hoodini/extra_tools/`
2. Follow interface: `run_*(input_files, output_dir, threads, **kwargs)`
3. Add CLI flag in `cli.py`
4. Update `write_data.py` to merge results

### Modifying pipeline stages

1. Each stage in `hoodini/pipeline/` or root module
2. Stages accept typed config, return DataFrames
3. Test with small fixtures before integration

### Working with configuration

1. Add new keys to `config/defaults.toml`
2. Update schema in `config/schema.py`
3. Document in `config/example.toml`
