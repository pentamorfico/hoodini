# Contributing to Hoodini

Thank you for your interest in contributing to Hoodini! This document covers our development standards and how to get started.

## Getting Started

### Development Setup

```bash
# Clone the repository
git clone https://github.com/pentamorfico/hoodini.git
cd hoodini

# Create development environment
mamba env create -f environment.yml
mamba activate hoodini

# Install in editable mode
pip install -e ".[dev]"
```

### Running Tests

```bash
pytest tests/
pytest tests/ -v --tb=short  # verbose with short tracebacks
```

## Code Style

We follow PEP 8 and the Google Python Style Guide with these specifics:

### Formatting

- **Line length:** 100 characters
- **Tools:** `black` + `isort` + `ruff`
- Run before committing: `black . && isort . && ruff check .`

### Imports

```python
# Order: stdlib, third-party, first-party (blank lines between)
import os
from pathlib import Path

import pandas as pd
from rich.console import Console

from hoodini.config import load_config
from hoodini.utils.core import validate_input
```

### Naming Conventions

| Type | Convention | Example |
|------|------------|---------|
| Modules | `snake_case` | `tree_builders.py` |
| Functions | `snake_case` verbs | `build_tree()`, `fetch_ipg()` |
| Classes | `PascalCase` nouns | `RuntimeConfig`, `AssemblyDownloader` |
| Constants | `UPPER_SNAKE_CASE` | `DEFAULT_THREADS` |
| Private | `_prefixed` | `_parse_header()` |

Reserve `run_*` for CLI-facing entry wrappers.

### Type Hints

```python
from __future__ import annotations
from pathlib import Path

def process_file(input_path: Path, threads: int = 4) -> pd.DataFrame:
    """Process input file and return results.
    
    Args:
        input_path: Path to input file.
        threads: Number of threads to use.
        
    Returns:
        DataFrame with processed results.
        
    Raises:
        FileNotFoundError: If input file doesn't exist.
    """
    ...
```

### Logging

```python
import logging

logger = logging.getLogger(__name__)

# Use appropriate levels
logger.debug("Processing file: %s", path)
logger.info("Completed %d assemblies", count)
logger.warning("Missing annotation for %s", gene_id)
logger.error("Failed to download: %s", url)
```

Avoid `print()` in library code—use logging instead. CLI UX messages are the exception.

### Error Handling

```python
# Good: specific exception with clear message
if not input_file.exists():
    raise FileNotFoundError(f"Input file not found: {input_file}")

# Bad: bare except
try:
    process()
except:  # Never do this
    pass
```

## Project Structure

```
hoodini/
├── cli.py           # Keep thin: parse args, build config, call pipeline
├── config/          # Configuration handling
├── pipeline/        # Core pipeline stages
├── utils/           # Focused helpers by concern
└── extra_tools/     # External tool wrappers
```

### Guidelines

- **I/O at the edges:** Core logic should be pure transformations
- **Typed configs:** Avoid ad-hoc dicts; use dataclasses
- **Small modules:** Split large files by concern
- **Test coverage:** Add tests for new functionality

## Pull Request Process

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/my-feature`
3. Make your changes with tests
4. Run linting: `black . && isort . && ruff check .`
5. Run tests: `pytest tests/`
6. Commit with clear message: `git commit -m "feat: add new annotation tool"`
7. Push and open a PR

### Commit Messages

Follow [Conventional Commits](https://www.conventionalcommits.org/):

```
feat: add PADLOC annotation support
fix: handle missing IPG records gracefully
docs: update installation instructions
refactor: split cli.py into modules
test: add fixtures for assembly parser
```

## Questions?

Open an issue or reach out to the maintainers.
