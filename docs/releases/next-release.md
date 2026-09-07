# Next release work log

## Scope and working agreement

- Working branch: `update/hoodini-next-release`.
- Starting commit: `c58abd7933b46822b29851872bded5b4b1e88549` (`main`).
- Started: 2026-09-07.
- Release number and publication date are intentionally unset.
- Keep implementation, validation results and remaining work in this file. Keep
  user-facing release notes in the root `CHANGELOG.md`.
- Commit related changes together and reference the original issue numbers.
- An implemented fix is not a closed issue or a published release. Completion
  requires its regression checks; release validation includes the full pipeline.

## Issue inventory

This inventory was checked against the 20 open issues and their comments. The
first implementation batch addresses data correctness and taxonomy. Subsequent
features remain candidates, not promised release content.

| Issue | Scope | Status / next action |
| --- | --- | --- |
| [#83](https://github.com/pentamorfico/hoodini/issues/83) | Multi-contig GFF selection | Implemented and regression-tested on this branch. Full-pipeline release check remains. |
| [#84](https://github.com/pentamorfico/hoodini/issues/84) | Mixed GBFF and GFF/FAA types | Implemented and regression-tested with real GFF/FAA and GBFF fixtures. |
| [#86](https://github.com/pentamorfico/hoodini/issues/86) | Heterogeneous contig Parquet schemas | Implemented and regression-tested across reader, updater and writer paths. Large-dataset validation remains. |
| [#81](https://github.com/pentamorfico/hoodini/issues/81) | NCBI domain taxonomy | Short-term compatibility fix implemented and regression-tested. Canonical-field migration and the assembly-summary parsing error in a comment remain separate work. |
| [#57](https://github.com/pentamorfico/hoodini/issues/57) | Optional database setup on later runs | Pending: inspect launcher and pipeline checks, including Colab. |
| [#75](https://github.com/pentamorfico/hoodini/issues/75) | DefenseFinder / CasFinder compatibility | Pending: verify installed tool and model compatibility before choosing constraints. |
| [#49](https://github.com/pentamorfico/hoodini/issues/49) | Run parameters in HTML | Pending: persist effective parameters and display them in the output. |
| [#78](https://github.com/pentamorfico/hoodini/issues/78) | Copy filtered table cells | Pending: reproduce in hoodini-viz. |
| [#76](https://github.com/pentamorfico/hoodini/issues/76) | Protein metadata links and copying | Pending: coordinate with the hoodini-viz table changes. |
| [#79](https://github.com/pentamorfico/hoodini/issues/79) | Search neighborhood metadata | Pending: implement and verify in hoodini-viz. |
| [#80](https://github.com/pentamorfico/hoodini/issues/80) | Mouse-wheel scrolling | Pending: implement and verify in hoodini-viz. |
| [#13](https://github.com/pentamorfico/hoodini/issues/13) | Presence/absence heatmap | Pending: define feature aggregation and viewer layout. |
| [#46](https://github.com/pentamorfico/hoodini/issues/46) | AAI/ANI trees and dereplication | Partially present on main: AAI/ANI trees; dereplication remains. |
| [#17](https://github.com/pentamorfico/hoodini/issues/17) | Cluster representatives by context completeness | Pending: specify deterministic selection and integrate with dereplication. |
| [#59](https://github.com/pentamorfico/hoodini/issues/59) | Streaming assemblies | Pending: define shared-assembly processing and interaction with --keep. |
| [#4](https://github.com/pentamorfico/hoodini/issues/4) | Resistance annotation | Pending: choose supported tool/database and acceptance fixtures. |
| [#3](https://github.com/pentamorfico/hoodini/issues/3) | BGC annotation | Pending: choose supported tool and acceptance fixtures. |
| [#64](https://github.com/pentamorfico/hoodini/issues/64) | Scheduled contig-table builds | Pending: establish build resources and publication destination. |
| [#85](https://github.com/pentamorfico/hoodini/issues/85) | Shared database updates | Pending: specify publication, trust and credentials. |
| [#77](https://github.com/pentamorfico/hoodini/issues/77) | Bioconda package | Pending: verify dependencies, recipe and external review. |

## Activity and evidence

### 2026-09-07: branch and baseline inspection

- Created the remote update branch from main and checked out the same base locally.
- Read `AGENTS.md`, `CONTRIBUTING.md`, the test configuration and the current CI.
- Confirmed the four priority defects are still present in the source.
- Found an additional writer hazard within #86: PyArrow infers dictionary columns
  from the first row, which can omit optional fields present in subsequent rows.
- Created this issue inventory and the Unreleased changelog.

### 2026-09-07: data-correctness batch

- #83: select a single GFF contig before computing coordinate or gene-count
  windows. Use the matching FNA record, infer the contig for unique protein IDs
  or single-contig inputs, and report ambiguous/missing contigs explicitly.
- #84: declare the nine GFF field types during CSV parsing, retaining score and
  phase as strings compatible with GBFF. Skip comment lines at parse time.
- #81: map domain into superkingdom only when the legacy rank has no value.
  Preserve existing output columns and do not treat root ranks as domains.
- #86: introduce a shared, lazy four-column contig scan. DuckDB combines schemas
  by name and supplies typed nulls for absent optional accession columns; the
  Polars fallback normalizes partitions individually before concatenation.
- Apply that scan to nuc2asmlen, direct pipeline lookup and missing-assembly
  detection. Keep the 4 GB DuckDB limit and close the touched query connections
  on errors as well as success. The updater now uses its supplied summary path.
- Preserve writer metadata from all buffered rows, including keys absent from
  the first row, and use stable Arrow types for the four core lookup fields.
- Added 30 regression cases in three test modules. The initial 28 cases produced
  22 failures and 6 passes against the unmodified implementation, reproducing
  the reported issues. Two additional selection guards were added afterward.
- Updated the user-facing Unreleased changelog with the implemented behavior.

### Execution environment and publication notes

- Local Git clone succeeded; terminal push lacked GitHub credentials. Repository
  writes therefore use the authenticated GitHub connector. The documentation
  baseline is commit `5781f344c0dc4b4831be9561ae581ff95c55bdcc`.
- Conda/Mamba was unavailable locally and its bootstrap download was interrupted
  by network approval. The full environment required by AGENTS.md for external
  bioinformatics tools has not been set up or claimed as validated.
- Used an isolated Python 3.12.13 environment for the Python regression suite,
  following the Python-only approach of the existing CI. No external analysis
  tools, live taxonomy downloads or complete NCBI databases were used by the new
  tests. NCBI lineage responses are deterministic fixtures.
- The initial orfipy build failed because the runtime requested missing clang;
  rebuilding with the installed GCC/G++ succeeded. An interrupted installer
  status was checked by an offline dependency check before running tests.
- Key validation versions: Polars 1.44.1, DuckDB 1.5.5, PyArrow 25.0.1,
  pytest 9.1.1, Biopython 1.88, gb-io 0.4.0, pyrodigal 3.7.1, orfipy 0.0.4,
  ETE3 3.1.3, AlphaFetcher 0.2.0, Black 24.10.0, isort 9.0.1, Ruff 0.16.6.

## Validation

| Check | Result |
| --- | --- |
| `pytest tests/unit tests/integration -o addopts='' -q --tb=short --disable-warnings` | 65 passed, 5 skipped. The five existing skips require contig or assembly-summary databases. |
| New regression cases | 30 passed (11 neighborhood, 7 taxonomy, 12 contig metadata). |
| Black and isort on the nine changed/new Python files | Passed. |
| Ruff on the four new Python files | Passed. |
| Ruff on the five modified Python files | Seven PLR0917 diagnostics, identical in count and rule to the original files at the starting commit. Pre-existing lint debt remains; the lint gate is not claimed to pass. |
| `git diff --check` | Passed. |
| Full Conda/Bioconda pipeline, live NCBI access, large-database performance | Not run; required before release. |
| Minimum supported Python/dependency versions | Not run locally; current-version results do not establish lower-bound compatibility. |

The GFF regression uses an actual GBFF written by Biopython and read by gb-io.
Contig tests exercise actual DuckDB/Polars queries and Parquet I/O, both partition
orders, optional columns absent from an entire dataset, forced Polars fallback,
pipeline enrichment and incremental writing. No synthetic test stubs replace
the parsers or database engines under test.

## Next implementation batch

1. Reproduce the assembly-summary parsing failure reported in #81's comment and
   add stable schema handling without silently publishing incomplete summaries.
2. Inspect #57 and #75 together using a real Conda environment and the Colab
   launcher, then test tool/model compatibility and optional database setup.
3. Add run provenance in the generated HTML (#49), excluding API keys.
4. Coordinate table and search work in hoodini-viz (#78, #76, #79, #80) and test
   the compiled viewer embedded in Hoodini output.

## Release checklist

- [x] Complete and regression-test the first data-correctness batch.
- [x] Review additional release scope against the issue inventory.
- [ ] Run the full Conda/Bioconda pipeline with representative local and downloaded data.
- [ ] Verify viewer changes in the generated standalone HTML, if included.
- [ ] Review the Unreleased notes and choose the release version.
- [ ] Merge, tag and publish the release when ready.
