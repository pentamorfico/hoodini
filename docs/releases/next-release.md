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
| [#83](https://github.com/pentamorfico/hoodini/issues/83) | Multi-contig GFF selection | In progress: isolate the requested contig before window extraction. |
| [#84](https://github.com/pentamorfico/hoodini/issues/84) | Mixed GBFF and GFF/FAA types | In progress: normalize GFF fields and test mixed inputs. |
| [#86](https://github.com/pentamorfico/hoodini/issues/86) | Heterogeneous contig Parquet schemas | In progress: normalize readers and writer, including absent optional columns. |
| [#81](https://github.com/pentamorfico/hoodini/issues/81) | NCBI domain taxonomy | In progress: accept domain while preserving legacy output. The assembly-summary parsing error in a comment needs separate evaluation. |
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

## Validation

No implementation validation has been run yet.

## Release checklist

- [ ] Complete and test the first data-correctness batch.
- [ ] Review additional release scope against the issue inventory.
- [ ] Run the full Conda/Bioconda pipeline with representative local and downloaded data.
- [ ] Verify viewer changes in the generated standalone HTML, if included.
- [ ] Review the Unreleased notes and choose the release version.
- [ ] Merge, tag and publish the release when ready.
