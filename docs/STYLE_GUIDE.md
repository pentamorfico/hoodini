# Hoodini Python Style Guide

This project follows PEP 8 and the Google Python Style Guide with the adjustments below.

## Imports
- Order: standard library, third-party, first-party (`hoodini`), with blank lines between groups.
- No wildcard imports; avoid unused imports; prefer module imports over `from module import *`.
- Lazy-import only to break cycles—document why when doing so.

## Naming
- Modules/packages: `snake_case` by domain (`tree_builders.py`, `nucleotide_links.py`).
- Functions/methods: `snake_case` verbs (`build_tree`, `fetch_ipg`); reserve `run_*` for CLI-facing entry wrappers.
- Classes: `PascalCase` nouns (`RuntimeConfig`, `AssemblyDownloader`).
- Constants: `UPPER_SNAKE_CASE`.
- Private helpers: prefix with `_` when not part of the public API.

## Types and Docstrings
- Enable `from __future__ import annotations` in new modules.
- All public functions/methods must be type hinted; prefer `Path` over raw strings for filesystem paths.
- Use Google-style docstrings with `Args`, `Returns`, and `Raises`; describe side effects clearly.

## Logging and Errors
- Library code uses `logging.getLogger(__name__)`; avoid `print` except for CLI UX.
- Raise specific exceptions; avoid bare `except` blocks. Fail fast with clear messages.

## Structure
- Keep I/O at the edges; core logic should be pure transformations.
- CLI files stay thin: parse args, build config, call orchestrated pipeline functions.
- Group helpers by concern in `hoodini/utils/` (e.g., `validation.py`, `polars_adapters.py`, `logging_utils.py`).

## Formatting
- Line length: 100; use `black` + `isort` (already configured). `ruff` is the primary linter.
- Comments are rare and explain intent, not mechanics.

## Testing
- Prefer pytest fixtures for file system or data setup.
- Add unit tests per pipeline stage plus a small end-to-end smoke test (guarded by an env flag for slow runs).
