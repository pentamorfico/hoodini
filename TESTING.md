# Hoodini Testing

## Structure

```
tests/
├── conftest.py              # Shared fixtures
├── unit/                    # pytest - fast, pure logic
│   ├── test_config.py
│   ├── test_windows.py
│   ├── test_cand_modes.py
│   └── test_helpers.py
├── integration/             # pytest - module wiring
│   └── test_runner.py
└── feature_coverage/        # Custom runner - real E2E
    ├── run_pipeline_matrix.py
    └── inputs/
```

## Running Tests

```bash
# Unit + integration (fast)
pytest tests/unit tests/integration -v

# With coverage
pytest tests/unit tests/integration --cov=hoodini

# Feature coverage matrix
python tests/feature_coverage/run_pipeline_matrix.py --all
python tests/feature_coverage/run_pipeline_matrix.py 1 2 3  # specific tests
python tests/feature_coverage/run_pipeline_matrix.py --coverage  # show matrix
python tests/feature_coverage/run_pipeline_matrix.py --download
python tests/feature_coverage/run_pipeline_matrix.py --utils
```

## Adding Tests

### Unit Test
```python
# tests/unit/test_<module>.py
def test_my_function():
    result = my_function(input)
    assert result == expected
```

### Integration Test
```python
# tests/integration/test_<feature>.py
def test_module_wiring(tmp_path):
    cfg = RuntimeConfig(...)
    # Test interaction between modules
```

### Feature Coverage Test
```python
# Add to TESTS in run_pipeline_matrix.py
TestCase(
    id=10,
    name="New Feature Test",
    covers=["new_flag", "another_flag"],
    config={...},
)
```

## CI

```yaml
# Fast (every commit)
pytest tests/unit tests/integration --cov=hoodini

# Full (nightly/release)
python tests/feature_coverage/run_pipeline_matrix.py --all
```
