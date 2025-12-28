# hoodini engineering & style guide

This guide collects conventions for organizing the codebase, naming things, structuring imports, and writing readable, well-typed Python. It draws primarily from [PEP 8](https://peps.python.org/pep-0008/) and the [Google Python Style Guide](https://google.github.io/styleguide/pyguide.html).

## Goals
- Make the repository easy to navigate by using predictable package boundaries and filenames.
- Encourage consistent naming (no surprise prefixes) across CLI entrypoints, pipelines, helpers, and tests.
- Prefer typed, documented functions with minimal side effects and explicit dependencies.
- Keep imports clean and deterministic to reduce merge noise.

## Repository layout recommendations
The current code is flat, with many single-purpose scripts. New or reorganized code should prefer a package layout with clear responsibilities:

```
hoodini/
  cli/                # click/typer/wrapper code that defines user-facing commands
  config/             # defaults, schema validation, config helpers
  pipelines/          # orchestration flows (formerly run_* or pipeline-like modules)
  io/                 # file readers/writers, download helpers, path utilities
  annotators/         # annotation-related modules (PHROG, padloc, deffinder, etc.)
  plotting/           # plotting and figure generation
  taxonomy/           # taxonomy helpers and data loaders
  utils/              # shared helpers (see below for structure)
  tests/              # mirrors package layout for unit/functional tests
```

When moving existing modules, prefer incremental changes (e.g., move `run_*` pipeline files into `pipelines/`, plotting helpers into `plotting/`, and download helpers into `io/`). Expose consolidated entrypoints through `hoodini/cli/__init__.py` and keep `__main__.py` thin.

### Utilities
Utility-style functions currently live in both `hoodini/utils.py` and `hoodini/utils/`. Consolidate helpers under a single `hoodini/utils/` package with focused modules:

- `utils/io.py`: filesystem, compression, and serialization helpers.
- `utils/paths.py`: path normalization and output directory creation.
- `utils/concurrency.py`: thread/process pool helpers, job submission wrappers.
- `utils/subprocess.py`: shared subprocess invocation logic.
- `utils/validation.py`: small validation functions used across the CLI and pipelines.

Avoid adding new top-level scripts for helpers; instead, expand these modules or create a clearly named sibling inside `utils/`.

## Naming conventions
- **Modules and packages**: `snake_case` filenames. Use package folders for logical groupings instead of long filenames.
- **Functions**: `snake_case`. Reserve the `run_*` prefix for high-level pipeline orchestration or CLI-invoked workflows. Internal helpers should use descriptive verbs (`load_config`, `build_tree`, `write_table`).
- **Classes**: `PascalCase`.
- **Constants**: `UPPER_SNAKE_CASE` defined near the top of a module.
- **Private/experimental APIs**: prefix with a single leading underscore. Avoid double-underscore name-mangling unless needed inside classes.
- **Variables**: short but descriptive names; prefer `path`, `output_dir`, `records`, `assemblies` over single-letter names (except simple loop indices).
- **CLI options**: mirror config keys; expose them through well-named functions inside `hoodini/cli/`.

## Typing
- Use type hints on all public functions, methods, and dataclass fields.
- Enable `from __future__ import annotations` in new modules to simplify forward references.
- Prefer `pathlib.Path` for filesystem paths and annotate accordingly.
- When returning structured data, use `TypedDict`/`NamedTuple`/`dataclass` instead of loose dictionaries.
- Avoid `Any` unless interacting with third-party libraries; document the reason in a short comment.

## Imports
- Order imports in three blocks separated by blank lines: standard library, third-party, then local package imports.
- Prefer absolute imports from `hoodini.` rather than relative imports between distant modules.
- Avoid wildcard imports. Import only what is used.
- Keep imports at the top of the file unless deferred imports materially reduce startup cost for the CLI.

## Docstrings and comments
- Use [Google-style docstrings](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings) for public functions, classes, and modules.
- Document arguments, return values, raised exceptions, and side effects (filesystem writes, network calls).
- Keep inline comments concise and explain *why* not *what* when the code is non-obvious.

## Error handling and logging
- Prefer `ValueError`, `KeyError`, `RuntimeError`, etc., over bare `Exception` when raising.
- Keep try/except blocks narrow and log context before re-raising.
- Route user-facing messages through `logging` instead of `print`, using module-level loggers named after `__name__`.

## CLI and orchestration patterns
- CLI functions should be thin: parse/validate inputs, construct configuration objects, and delegate to pipeline functions (e.g., `pipelines/run_pipeline.py:run_pipeline`).
- Pipeline functions should accept explicit arguments/config objects rather than reading globals or mutating module state.
- Prefer returning results instead of writing directly to disk inside deep utilities; let orchestration layers decide persistence.

## Tests
- Mirror the package layout under `tests/` (e.g., `tests/pipelines/test_run_pipeline.py`).
- Keep fixtures in `tests/fixtures/` or `tests/resources/` and avoid cross-test side effects.
- Use descriptive test names (`test_build_tree_rejects_empty_input`) and arrange tests in Arrange/Act/Assert order with clear comments if needed.

## Tooling and CI
- Run `ruff check`, `black`, and `isort` locally before opening a PR. The CI workflow enforces these tools **only on Python files changed in the PR** to prevent legacy formatting issues from blocking documentation-only or small fixes.
- Use `mypy` in lenient mode (current CI behavior) and tighten settings as modules gain full typing coverage.
- Add or update tests whenever changing behavior; CI runs `pytest` when the `tests/` folder is present.

## Incremental cleanup roadmap
1. Migrate shared helpers into the structured `hoodini/utils/` package and update imports.
2. Group CLI code under `hoodini/cli/` and delegate to pipeline modules in `hoodini/pipelines/`.
3. Standardize function prefixes: `run_*` only for top-level workflows; helper verbs elsewhere.
4. Add type hints and Google-style docstrings to public modules as they are touched.
5. Introduce linting/formatting tooling (e.g., `ruff`, `black`, `mypy`) and gradually expand CI enforcement from "changed files only" to full-tree checks as the backlog is addressed.

Applying these conventions incrementally will give the codebase a predictable shape while minimizing disruption to existing workflows.
