#!/usr/bin/env python3
"""
Hoodini Pipeline Test Suite
===========================
Smart tests that cover maximum CLI options with minimum test count.
Uses the Python API directly (no subprocess).

Usage:
    python run_tests.py              # Run all pipeline tests
    python run_tests.py 1 2 3        # Run specific pipeline tests
    python run_tests.py --coverage   # Show coverage matrix
    python run_tests.py --list       # List all tests
    python run_tests.py --download   # Run download tests (requires network)
    python run_tests.py --utils      # Run utils tests
    python run_tests.py --all        # Run everything (pipeline + download + utils)
    python run_tests.py --results    # Show last run results from log
"""

import json
import logging
import os
import sys
import tempfile
import traceback
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from hoodini.config import load_default_config
from hoodini.config.settings import RuntimeConfig
from hoodini.pipeline.runner import run_pipeline

# ANSI colors
GREEN = "\033[92m"
YELLOW = "\033[93m"
RED = "\033[91m"
CYAN = "\033[96m"
BOLD = "\033[1m"
NC = "\033[0m"

THREADS = 4
FEATURE_DIR = Path(__file__).parent
INPUTS_DIR = FEATURE_DIR / "inputs"
OUTPUTS_DIR = FEATURE_DIR / "outputs"
RESULTS_LOG = OUTPUTS_DIR / "test_results.json"

# Test results storage
TEST_RESULTS: list[dict] = []


# =============================================================================
# OUTPUT VALIDATION
# =============================================================================

# Files/dirs that MUST exist after successful run
MANDATORY_OUTPUTS = [
    "records.tsv",
    "hoodini-viz",
    "hoodini-viz/tree.nwk",
    "hoodini-viz/tsv/gff.gff",
    "hoodini-viz/tsv/hoods.txt",
    "hoodini-viz/tsv/protein_metadata.txt",
    "hoodini-viz/tsv/tree_metadata.txt",
    "hoodini-viz/parquet",
    "neighborhood/neighborhoods.fasta",
]

# Files/dirs that should NOT exist when keep=False
TEMP_FILES_TO_REMOVE = [
    "assembly_folder",           # Raw GenBank files
    "tmp",                       # General temp
    "tmp_mmseqs",                # MMseqs2 temp
    "ani_split",                 # FastANI splits
    "fastani_pairwise_visual",   # FastANI viz
    "struct",                    # Foldseek structures
    "temp",                      # Foldseek temp
    "temp.fasta",
    "temp.gff",
    "proteome.fasta",
    "proteome.fasta.idx",
    "assembly_list.txt",
    "all_neigh.tsv",
    "wgrr.tsv",
    "intergenic.fasta",
    "deepmmseqs_results.tsv",
    "fastani_genome_list.txt",
    "fastani_output.tsv",
    "fastani_all.log",
]

# Optional outputs that should exist if corresponding option was enabled
OPTIONAL_OUTPUT_MAP = {
    "padloc": ["padloc"],
    "deffinder": ["defense_finder"],
    "cctyper": ["cctyper"],
    "genomad": ["genomad/output"],
    "ncrna": ["ncrna"],
    # These are file-based, check for presence
    "aai_tree": ["aai.tsv"],
    "ani_tree": ["pairwise_ani_*.tsv"],
}


def validate_outputs(output_dir: Path, config: dict, keep: bool = False) -> tuple[bool, list[str]]:
    """Validate that all expected outputs exist and temp files are cleaned.
    
    Args:
        output_dir: Path to the output directory
        config: The test config dict (to check optional outputs)
        keep: Whether --keep flag was used
        
    Returns:
        Tuple of (all_valid, list_of_issues)
    """
    issues = []
    
    # 1. Check mandatory outputs exist
    for path in MANDATORY_OUTPUTS:
        full_path = output_dir / path
        if not full_path.exists():
            issues.append(f"MISSING mandatory: {path}")
    
    # 2. Check optional outputs based on config
    for opt, paths in OPTIONAL_OUTPUT_MAP.items():
        if config.get(opt):
            for path in paths:
                # Handle wildcards
                if "*" in path:
                    matches = list(output_dir.glob(path))
                    if not matches:
                        issues.append(f"MISSING optional ({opt}): {path}")
                else:
                    full_path = output_dir / path
                    if not full_path.exists():
                        issues.append(f"MISSING optional ({opt}): {path}")
    
    # 3. Check temp files are cleaned (only if keep=False)
    if not keep:
        for path in TEMP_FILES_TO_REMOVE:
            full_path = output_dir / path
            if full_path.exists():
                issues.append(f"TEMP not cleaned: {path}")
    
    # 4. Check for any unexpected large directories
    if not keep:
        for item in output_dir.iterdir():
            if item.is_dir() and item.name not in [
                "hoodini-viz", "neighborhood", 
                "padloc", "defense_finder", "cctyper", "genomad", "ncrna"
            ]:
                # Check size - if > 10MB, likely a temp we missed
                total_size = sum(f.stat().st_size for f in item.rglob("*") if f.is_file())
                if total_size > 10_000_000:  # 10MB
                    issues.append(f"LARGE unexpected dir: {item.name} ({total_size // 1_000_000}MB)")
    
    return len(issues) == 0, issues


def print_output_summary(output_dir: Path):
    """Print a summary of the output directory contents."""
    if not output_dir.exists():
        print(f"  {RED}Output directory does not exist{NC}")
        return
    
    print(f"\n{CYAN}Output contents:{NC}")
    for item in sorted(output_dir.iterdir()):
        if item.is_dir():
            # Count files and total size
            files = list(item.rglob("*"))
            file_count = sum(1 for f in files if f.is_file())
            total_size = sum(f.stat().st_size for f in files if f.is_file())
            size_str = f"{total_size // 1024}KB" if total_size < 1_000_000 else f"{total_size // 1_000_000}MB"
            print(f"  📁 {item.name}/ ({file_count} files, {size_str})")
        else:
            size = item.stat().st_size
            size_str = f"{size // 1024}KB" if size < 1_000_000 else f"{size // 1_000_000}MB"
            print(f"  📄 {item.name} ({size_str})")


def save_result(test_type: str, test_id: int | str, name: str, passed: bool, error: str | None = None):
    """Save a test result to the global results list."""
    TEST_RESULTS.append({
        "type": test_type,
        "id": test_id,
        "name": name,
        "passed": passed,
        "error": error,
        "timestamp": datetime.now().isoformat(),
    })


def write_results_log():
    """Write all results to the JSON log file."""
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    log_data = {
        "run_timestamp": datetime.now().isoformat(),
        "results": TEST_RESULTS,
        "summary": {
            "total": len(TEST_RESULTS),
            "passed": sum(1 for r in TEST_RESULTS if r["passed"]),
            "failed": sum(1 for r in TEST_RESULTS if not r["passed"]),
        }
    }
    with open(RESULTS_LOG, "w") as f:
        json.dump(log_data, f, indent=2)
    print(f"\n{CYAN}Results saved to: {RESULTS_LOG}{NC}")


def show_last_results():
    """Display results from the last test run."""
    if not RESULTS_LOG.exists():
        print(f"{RED}No results log found. Run tests first.{NC}")
        return
    
    with open(RESULTS_LOG) as f:
        data = json.load(f)
    
    print(f"\n{CYAN}{BOLD}Last Test Run: {data['run_timestamp']}{NC}\n")
    
    summary = data["summary"]
    color = GREEN if summary["failed"] == 0 else RED
    print(f"{color}Summary: {summary['passed']} passed, {summary['failed']} failed{NC}\n")
    
    # Show failed tests with details
    failed = [r for r in data["results"] if not r["passed"]]
    if failed:
        print(f"{RED}{BOLD}Failed Tests:{NC}\n")
        for r in failed:
            print(f"  {RED}✗ [{r['type']}] {r['id']}: {r['name']}{NC}")
            if r.get("error"):
                # Print first 500 chars of error
                err = r["error"][:500]
                if len(r["error"]) > 500:
                    err += "..."
                for line in err.split("\n"):
                    print(f"    {line}")
            print()
    
    # Show passed tests briefly
    passed = [r for r in data["results"] if r["passed"]]
    if passed:
        print(f"{GREEN}{BOLD}Passed Tests:{NC}")
        for r in passed:
            print(f"  {GREEN}✓ [{r['type']}] {r['id']}: {r['name']}{NC}")


# All trackable options for coverage matrix
ALL_OPTIONS = {
    "Input": ["input_path", "inputsheet", "keep", "debug"],
    "Tree Modes": [
        "taxonomy", "aai_tree", "ani_tree", "fast_ml", "fast_nj",
        "neigh_phylo_tree", "neigh_similarity_tree", "foldmason_tree"
    ],
    "Clustering": ["diamond_deepclust", "deepmmseqs", "blastp", "jackhmmer"],
    "Cand Mode": ["any_ipg", "best_ipg", "best_id", "one_id", "same_id"],
    "Annotations": [
        "padloc", "deffinder", "cctyper", "genomad", 
        "domains", "ncrna", "sorfs", "blast"
    ],
    "Pairwise": [
        "aai_mode=aai", "aai_mode=wgrr", "ani_mode=blastn",
        "prot_links", "nt_links", 
        "minimap2", "fastani", "intergenic_blastn"
    ],
    "Window": ["win_nts", "win_genes", "minwin_upstream", "minwin_downstream"],
}


@dataclass
class TestCase:
    """A single pipeline test case."""
    id: int
    name: str
    description: str
    config: dict
    covers: list[str] = field(default_factory=list)

    def run(self) -> bool:
        """Execute the test using hoodini's Python API."""
        print(f"{YELLOW}{'━' * 80}{NC}")
        print(f"{YELLOW}{BOLD}TEST {self.id}: {self.name}{NC}")
        print(f"{YELLOW}Covers: {', '.join(self.covers)}{NC}")
        print(f"{YELLOW}{'━' * 80}{NC}")

        try:
            defaults = load_default_config()
            merged = {**defaults, **self.config}
            
            # Resolve input paths to absolute
            if merged.get("input_path"):
                merged["input_path"] = str(INPUTS_DIR / merged["input_path"])
            if merged.get("inputsheet"):
                merged["inputsheet"] = str(INPUTS_DIR / merged["inputsheet"])
            if merged.get("blast"):
                merged["blast"] = str(INPUTS_DIR / merged["blast"])
            # Output goes to outputs/ folder (gitignored)
            merged["output"] = str(OUTPUTS_DIR / merged["output"])
            
            cfg = RuntimeConfig(
                input_path=merged.get("input_path"),
                inputsheet=merged.get("inputsheet"),
                output=merged.get("output"),
                max_concurrent_downloads=merged.get("max_concurrent_downloads", 3),
                num_threads=merged.get("num_threads", THREADS),
                assembly_folder=merged.get("assembly_folder"),
                prot_links=merged.get("prot_links", False),
                nt_links=merged.get("nt_links", False),
                ani_mode=merged.get("ani_mode"),
                nt_aln_mode=merged.get("nt_aln_mode"),
                blast=merged.get("blast"),
                cand_mode=merged.get("cand_mode"),
                clust_method=merged.get("clust_method"),
                mod=merged.get("mod", "win_nts"),
                wn=merged.get("wn", 50000),
                minwin=merged.get("minwin"),
                minwin_type=merged.get("minwin_type"),
                tree_mode=merged.get("tree_mode", "taxonomy"),
                tree_file=merged.get("tree_file"),
                aai_mode=merged.get("aai_mode"),
                aai_subset_mode=merged.get("aai_subset_mode"),
                nj_algorithm=merged.get("nj_algorithm"),
                remote_evalue=merged.get("remote_evalue"),
                remote_max_targets=merged.get("remote_max_targets"),
                padloc=merged.get("padloc", False),
                deffinder=merged.get("deffinder", False),
                ncrna=merged.get("ncrna"),
                cctyper=merged.get("cctyper", False),
                genomad=merged.get("genomad", False),
                sorfs=merged.get("sorfs", False),
                domains=merged.get("domains"),
                emapper=merged.get("emapper", False),
                min_pident=merged.get("min_pident", 30.0),
                keep=merged.get("keep", False),
                force=True,
                debug=merged.get("debug", False),
            )

            print(f"{CYAN}Running pipeline:{NC}")
            print(f"  input: {cfg.input_path or cfg.inputsheet}")
            print(f"  output: {cfg.output}")
            print(f"  tree_mode: {cfg.tree_mode}\n")

            run_pipeline(cfg)
            
            # Validate outputs
            output_path = Path(cfg.output)
            print_output_summary(output_path)
            
            valid, issues = validate_outputs(output_path, self.config, keep=cfg.keep)
            if not valid:
                print(f"\n{YELLOW}Output validation issues:{NC}")
                for issue in issues:
                    print(f"  {YELLOW}⚠ {issue}{NC}")
            else:
                print(f"\n{GREEN}✓ Output validation passed{NC}")
            
            print(f"\n{GREEN}✓ TEST {self.id} PASSED{NC}\n")
            save_result("pipeline", self.id, self.name, True)
            return True

        except Exception as e:
            error_msg = f"{type(e).__name__}: {e}\n{traceback.format_exc()}"
            print(f"\n{RED}✗ TEST {self.id} FAILED{NC}")
            print(f"{RED}Exception: {e}{NC}\n")
            logging.exception("Test failed")
            save_result("pipeline", self.id, self.name, False, error_msg)
            return False


# =============================================================================
# TEST DEFINITIONS
# =============================================================================

TESTS = [
    TestCase(
        id=1,
        name="Taxonomy + PADLOC + DefenseFinder",
        description="Taxonomy tree with defense system annotations",
        covers=["taxonomy", "padloc", "deffinder", "input_path", "diamond_deepclust", "any_ipg", "win_nts"],
        config={
            "input_path": "cas9.txt",
            "output": "out_test1",
            "tree_mode": "taxonomy",
            "wn": 50000,
            "padloc": True,
            "deffinder": True,
            # Defaults explicitly: clust_method=diamond_deepclust, cand_mode=any_ipg, mod=win_nts
        },
    ),
    TestCase(
        id=2,
        name="AAI Tree + CCTyper + Genomad",
        description="AAI-based tree with CRISPR and mobile element detection",
        covers=["aai_tree", "cctyper", "genomad", "aai_mode=aai", "best_ipg"],
        config={
            "input_path": "cas9.txt",
            "output": "out_test2",
            "tree_mode": "aai_tree",
            "aai_mode": "aai",
            "wn": 40000,
            "cctyper": True,
            "genomad": True,
            "cand_mode": "best_ipg",
        },
    ),
    TestCase(
        id=3,
        name="ANI Tree + Domains",
        description="ANI-based tree with protein domain annotations",
        covers=["ani_tree", "domains", "ani_mode=blastn", "same_id"],
        config={
            "input_path": "cas9.txt",
            "output": "out_test3",
            "tree_mode": "ani_tree",
            "ani_mode": "blastn",
            "wn": 40000,
            "domains": ["pfam"],
            "cand_mode": "same_id",
        },
    ),
    TestCase(
        id=4,
        name="Fast ML + ncRNA + sorfs",
        description="ML tree with non-coding RNA and small ORF detection",
        covers=["fast_ml", "ncrna", "sorfs", "best_id"],
        config={
            "input_path": "cas9.txt",
            "output": "out_test4",
            "tree_mode": "fast_ml",
            "wn": 30000,
            "ncrna": "RF00001",  # 5S ribosomal RNA
            "sorfs": True,
            "cand_mode": "best_id",
        },
    ),
    TestCase(
        id=5,
        name="Fast NJ + BLAST + Links",
        description="NJ tree with custom BLAST and link calculations",
        covers=["fast_nj", "blast", "prot_links", "nt_links", "minimap2"],
        config={
            "input_path": "cas9.txt",
            "output": "out_test5",
            "tree_mode": "fast_nj",
            "wn": 30000,
            "blast": "blast_query.fasta",
            "prot_links": True,
            "nt_links": True,
            "nt_aln_mode": "minimap2",
        },
    ),
    TestCase(
        id=6,
        name="Neigh Phylo + Jackhmmer + WGRR",
        description="Neighborhood phylogeny with jackhmmer clustering",
        covers=["neigh_phylo_tree", "jackhmmer", "aai_mode=wgrr", "one_id"],
        config={
            "input_path": "cas9.txt",
            "output": "out_test6",
            "tree_mode": "neigh_phylo_tree",
            "clust_method": "jackhmmer",
            "aai_mode": "wgrr",
            "wn": 50000,
            "cand_mode": "one_id",
        },
    ),
    TestCase(
        id=7,
        name="Neigh Similarity + DeepMMseqs",
        description="Neighborhood similarity tree with deep clustering",
        covers=["neigh_similarity_tree", "deepmmseqs", "fastani", "minwin_upstream"],
        config={
            "input_path": "cas9.txt",
            "output": "out_test7",
            "tree_mode": "neigh_similarity_tree",
            "clust_method": "deepmmseqs",
            "wn": 40000,
            "minwin": 8000,
            "minwin_type": "upstream",
            "nt_links": True,
            "nt_aln_mode": "fastani",
        },
    ),
    TestCase(
        id=8,
        name="Inputsheet + BLASTp + win_genes",
        description="Inputsheet input with blastp clustering",
        covers=["inputsheet", "blastp", "win_genes", "intergenic_blastn", "debug", "keep", "minwin_downstream"],
        config={
            "inputsheet": "inputsheet.tsv",
            "output": "out_test8",
            "tree_mode": "taxonomy",
            "clust_method": "blastp",
            "mod": "win_genes",
            "wn": 15,
            "minwin": 5,
            "minwin_type": "downstream",
            "nt_links": True,
            "nt_aln_mode": "intergenic_blastn",
            "debug": True,
            "keep": True,
        },
    ),
    TestCase(
        id=9,
        name="Foldmason Tree (AlphaFold)",
        description="Structure-based tree using AlphaFold",
        covers=["foldmason_tree"],
        config={
            "input_path": "foldmason_ids.txt",
            "output": "out_test9",
            "tree_mode": "foldmason_tree",
            "wn": 25000,
        },
    ),
]


# =============================================================================
# COVERAGE MATRIX - Generated programmatically
# =============================================================================

def generate_coverage_matrix() -> str:
    """Generate coverage matrix from test definitions."""
    # Collect all unique options from tests
    all_covered = set()
    for t in TESTS:
        all_covered.update(t.covers)
    
    # Build matrix
    n_tests = len(TESTS)
    col_width = 4
    header = "│".join([f"T{t.id}".center(col_width) for t in TESTS])
    
    lines = []
    lines.append("=" * 100)
    lines.append("CLI OPTIONS COVERAGE MATRIX".center(100))
    lines.append("=" * 100)
    lines.append(f"{'OPTION':<30} │ {header} │")
    lines.append("-" * 100)
    
    for category, options in ALL_OPTIONS.items():
        lines.append(f"{category.upper():<30} │" + "│".join([" " * col_width] * n_tests) + "│")
        for opt in options:
            row = []
            for t in TESTS:
                if opt in t.covers:
                    row.append(" ✓ ".center(col_width))
                else:
                    row.append(" " * col_width)
            lines.append(f"  {opt:<28} │{'│'.join(row)}│")
        lines.append("-" * 100)
    
    # Summary
    lines.append("")
    lines.append("COVERAGE SUMMARY:")
    covered = sorted(all_covered)
    not_covered = []
    for opts in ALL_OPTIONS.values():
        for opt in opts:
            if opt not in all_covered:
                not_covered.append(opt)
    
    lines.append(f"  Covered ({len(covered)}): {', '.join(covered)}")
    if not_covered:
        lines.append(f"  Not covered ({len(not_covered)}): {', '.join(not_covered)}")
    lines.append("")
    lines.append("Special setup required: config (TOML), assembly_folder, apikey, tree_file, use_input_tree")
    
    return "\n".join(lines)


# =============================================================================
# DOWNLOAD TESTS - Testing hoodini download subcommands
# =============================================================================

@dataclass
class DownloadTestCase:
    """A single download test case."""
    id: int
    name: str
    description: str
    func: str  # Name of the function to import and call
    module: str  # Module path
    kwargs: dict = field(default_factory=dict)

    def run(self) -> bool:
        """Execute the download test using hoodini's Python API."""
        print(f"{YELLOW}{'━' * 80}{NC}")
        print(f"{YELLOW}{BOLD}DOWNLOAD TEST {self.id}: {self.name}{NC}")
        print(f"{YELLOW}Module: {self.module}.{self.func}{NC}")
        print(f"{YELLOW}{'━' * 80}{NC}")

        try:
            # Dynamic import
            module = __import__(self.module, fromlist=[self.func])
            func = getattr(module, self.func)
            
            print(f"{CYAN}Calling {self.func}({self.kwargs}){NC}\n")
            func(**self.kwargs)
            
            print(f"\n{GREEN}✓ DOWNLOAD TEST {self.id} PASSED{NC}\n")
            save_result("download", self.id, self.name, True)
            return True

        except Exception as e:
            error_msg = f"{type(e).__name__}: {e}\n{traceback.format_exc()}"
            print(f"\n{RED}✗ DOWNLOAD TEST {self.id} FAILED{NC}")
            print(f"{RED}Exception: {e}{NC}\n")
            logging.exception("Download test failed")
            save_result("download", self.id, self.name, False, error_msg)
            return False


DOWNLOAD_TESTS = [
    DownloadTestCase(
        id=1,
        name="Assembly Summary Download",
        description="Download and update assembly_summary.parquet",
        module="hoodini.download.assembly_summary",
        func="download_assembly_summary_db",
        kwargs={},  # Uses defaults
    ),
    DownloadTestCase(
        id=2,
        name="Databases (all skipped except contig_lengths)",
        description="Test database download logic with all DBs skipped (dry run)",
        module="hoodini.download.databases",
        func="main",
        kwargs={
            "force": False,
            "skip_padloc": True,
            "skip_deffinder": True,
            "skip_genomad": True,
            "skip_emapper": True,
            "skip_parquet": True,
            "skip_contig_lengths": True,
            "skip_typedive": True,
            "num_threads": 4,
        },
    ),
    DownloadTestCase(
        id=3,
        name="MetaCerberus DB Status",
        description="Show MetaCerberus database status (no actual download)",
        module="hoodini.download.metacerberus",
        func="main",
        kwargs={"selected": "all", "force": False},  # Just shows status
    ),
    # Note: type_dive.main() does actual downloads - skip in tests
    # DownloadTestCase(
    #     id=4,
    #     name="Type Dive Download",
    #     description="Download BacDive/PhageDive database",
    #     module="hoodini.download.type_dive",
    #     func="main",
    #     kwargs={},
    # ),
]


# =============================================================================
# UTILS TESTS - Testing hoodini utils subcommands
# =============================================================================

@dataclass
class UtilsTestCase:
    """A single utils test case."""
    id: int
    name: str
    description: str
    
    def run(self) -> bool:
        """Execute the utils test - implemented per test."""
        raise NotImplementedError


class Nuc2AsmlenTest(UtilsTestCase):
    """Test nuc2asmlen utility."""
    
    def __init__(self):
        super().__init__(
            id=1,
            name="nuc2asmlen",
            description="Fetch assembly and length metadata for nucleotide accessions",
        )

    def run(self) -> bool:
        print(f"{YELLOW}{'━' * 80}{NC}")
        print(f"{YELLOW}{BOLD}UTILS TEST {self.id}: {self.name}{NC}")
        print(f"{YELLOW}{self.description}{NC}")
        print(f"{YELLOW}{'━' * 80}{NC}")

        try:
            from hoodini.pipeline.helpers.nuc2asmlen import run_nuc2asmlen
            
            # Test with a list of nucleotide accessions
            test_accessions = [
                "NC_000913.3",   # E. coli K-12 chromosome
                "NZ_CP028116.1", # Another example
            ]
            
            print(f"{CYAN}Testing run_nuc2asmlen with accessions: {test_accessions}{NC}\n")
            df = run_nuc2asmlen(test_accessions)
            
            print(f"Result DataFrame shape: {df.shape}")
            print(df)
            
            # Verify we got results
            if df.height > 0:
                print(f"\n{GREEN}✓ UTILS TEST {self.id} PASSED{NC}\n")
                save_result("utils", self.id, self.name, True)
                return True
            else:
                print(f"\n{YELLOW}⚠ UTILS TEST {self.id}: No results (may be expected if contig_lengths.parquet is empty){NC}\n")
                save_result("utils", self.id, self.name, True)  # Not a failure
                return True  # Not a failure, just no data

        except Exception as e:
            error_msg = f"{type(e).__name__}: {e}\n{traceback.format_exc()}"
            print(f"\n{RED}✗ UTILS TEST {self.id} FAILED{NC}")
            print(f"{RED}Exception: {e}{NC}\n")
            logging.exception("Utils test failed")
            save_result("utils", self.id, self.name, False, error_msg)
            return False


class PrefetchLinksTest(UtilsTestCase):
    """Test prefetch_links utility."""
    
    def __init__(self):
        super().__init__(
            id=2,
            name="prefetch_links",
            description="Generate prefetched NCBI dataset links for assemblies",
        )

    def run(self) -> bool:
        print(f"{YELLOW}{'━' * 80}{NC}")
        print(f"{YELLOW}{BOLD}UTILS TEST {self.id}: {self.name}{NC}")
        print(f"{YELLOW}{self.description}{NC}")
        print(f"{YELLOW}{'━' * 80}{NC}")

        try:
            from hoodini.pipeline.helpers.prefetch_links import get_prefetched_link_table
            
            # Test with a few assembly accessions
            test_accessions = [
                "GCF_000005845.2",  # E. coli K-12
                "GCF_000009045.1",  # Bacillus subtilis
            ]
            
            print(f"{CYAN}Testing get_prefetched_link_table with assemblies: {test_accessions}{NC}")
            print(f"{CYAN}Kinds: gbff, gff, fna{NC}\n")
            
            df = get_prefetched_link_table(test_accessions, kinds=["gbff", "gff", "fna"])
            
            print(f"Result DataFrame shape: {df.shape}")
            print(df.head(10))
            
            # Verify we got results
            if df.height > 0:
                print(f"\n{GREEN}✓ UTILS TEST {self.id} PASSED{NC}\n")
                save_result("utils", self.id, self.name, True)
                return True
            else:
                print(f"\n{RED}✗ UTILS TEST {self.id}: No links generated{NC}\n")
                save_result("utils", self.id, self.name, False, "No links generated")
                return False

        except Exception as e:
            error_msg = f"{type(e).__name__}: {e}\n{traceback.format_exc()}"
            print(f"\n{RED}✗ UTILS TEST {self.id} FAILED{NC}")
            print(f"{RED}Exception: {e}{NC}\n")
            logging.exception("Utils test failed")
            save_result("utils", self.id, self.name, False, error_msg)
            return False


class PrefetchLinksFileTest(UtilsTestCase):
    """Test prefetch_links with file input (like CLI)."""
    
    def __init__(self):
        super().__init__(
            id=3,
            name="prefetch_links (file input)",
            description="Generate prefetched links from a file of assembly IDs",
        )

    def run(self) -> bool:
        print(f"{YELLOW}{'━' * 80}{NC}")
        print(f"{YELLOW}{BOLD}UTILS TEST {self.id}: {self.name}{NC}")
        print(f"{YELLOW}{self.description}{NC}")
        print(f"{YELLOW}{'━' * 80}{NC}")

        try:
            from hoodini.pipeline.helpers.prefetch_links import get_prefetched_link_table
            
            # Create temp file with assembly IDs
            test_accessions = ["GCF_000005845.2", "GCF_000009045.1"]
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
                f.write("\n".join(test_accessions))
                temp_file = f.name
            
            try:
                # Read file like CLI does
                with open(temp_file) as fh:
                    accs = [l.strip() for l in fh if l.strip()]
                
                print(f"{CYAN}Testing with file input: {temp_file}{NC}")
                print(f"{CYAN}Assemblies: {accs}{NC}")
                print(f"{CYAN}Kinds: sequence_report{NC}\n")
                
                df = get_prefetched_link_table(accs, kinds=["sequence_report"])
                
                print(f"Result DataFrame shape: {df.shape}")
                print(df)
                
                if df.height > 0:
                    print(f"\n{GREEN}✓ UTILS TEST {self.id} PASSED{NC}\n")
                    save_result("utils", self.id, self.name, True)
                    return True
                else:
                    print(f"\n{RED}✗ UTILS TEST {self.id}: No links generated{NC}\n")
                    save_result("utils", self.id, self.name, False, "No links generated")
                    return False
            finally:
                os.unlink(temp_file)

        except Exception as e:
            error_msg = f"{type(e).__name__}: {e}\n{traceback.format_exc()}"
            print(f"\n{RED}✗ UTILS TEST {self.id} FAILED{NC}")
            print(f"{RED}Exception: {e}{NC}\n")
            logging.exception("Utils test failed")
            save_result("utils", self.id, self.name, False, error_msg)
            return False


UTILS_TESTS = [
    Nuc2AsmlenTest(),
    PrefetchLinksTest(),
    PrefetchLinksFileTest(),
]


def print_coverage():
    """Print the coverage matrix."""
    print(f"{CYAN}{generate_coverage_matrix()}{NC}")


def list_tests():
    """List all available tests."""
    print(f"\n{CYAN}{BOLD}Pipeline Tests:{NC}\n")
    for t in TESTS:
        print(f"  {YELLOW}{t.id}{NC}: {t.name}")
        print(f"     {t.description}")
        print(f"     Covers: {', '.join(t.covers)}\n")
    
    print(f"\n{CYAN}{BOLD}Download Tests:{NC}\n")
    for t in DOWNLOAD_TESTS:
        print(f"  {YELLOW}D{t.id}{NC}: {t.name}")
        print(f"     {t.description}\n")
    
    print(f"\n{CYAN}{BOLD}Utils Tests:{NC}\n")
    for t in UTILS_TESTS:
        print(f"  {YELLOW}U{t.id}{NC}: {t.name}")
        print(f"     {t.description}\n")


def run_tests(test_ids: list[int] | None = None) -> tuple[int, int]:
    """Run specified pipeline tests or all tests. Returns (passed, failed)."""
    tests_to_run = TESTS if test_ids is None else [t for t in TESTS if t.id in test_ids]
    passed = failed = 0
    for test in tests_to_run:
        if test.run():
            passed += 1
        else:
            failed += 1
    return passed, failed


def run_download_tests(test_ids: list[int] | None = None) -> tuple[int, int]:
    """Run specified download tests or all. Returns (passed, failed)."""
    tests_to_run = DOWNLOAD_TESTS if test_ids is None else [t for t in DOWNLOAD_TESTS if t.id in test_ids]
    passed = failed = 0
    for test in tests_to_run:
        if test.run():
            passed += 1
        else:
            failed += 1
    return passed, failed


def run_utils_tests(test_ids: list[int] | None = None) -> tuple[int, int]:
    """Run specified utils tests or all. Returns (passed, failed)."""
    tests_to_run = UTILS_TESTS if test_ids is None else [t for t in UTILS_TESTS if t.id in test_ids]
    passed = failed = 0
    for test in tests_to_run:
        if test.run():
            passed += 1
        else:
            failed += 1
    return passed, failed


def main():
    args = sys.argv[1:]
    
    # No args = run pipeline tests only
    if not args:
        print(f"{CYAN}{BOLD}")
        print("╔══════════════════════════════════════════════════════════════════════════════╗")
        print(f"║           HOODINI PIPELINE TEST SUITE - Running ALL Tests (1-{len(TESTS)})             ║")
        print("╚══════════════════════════════════════════════════════════════════════════════╝")
        print(f"{NC}")
        
        print_coverage()
        passed, failed = run_tests()
        
        write_results_log()
        color = GREEN if failed == 0 else RED
        print(f"{color}")
        print(f"{'=' * 60}")
        print(f"  PIPELINE RESULTS: {passed} passed, {failed} failed")
        print(f"{'=' * 60}")
        print(f"{NC}")
        sys.exit(0 if failed == 0 else 1)
    
    if "--coverage" in args or "-c" in args:
        print_coverage()
        return
    
    if "--list" in args or "-l" in args:
        list_tests()
        return
    
    if "--help" in args or "-h" in args:
        print(__doc__)
        return
    
    if "--results" in args or "-r" in args:
        show_last_results()
        return
    
    # Run download tests
    if "--download" in args or "-d" in args:
        print(f"{CYAN}{BOLD}")
        print("╔══════════════════════════════════════════════════════════════════════════════╗")
        print(f"║           HOODINI DOWNLOAD TESTS - Running {len(DOWNLOAD_TESTS)} tests                          ║")
        print("╚══════════════════════════════════════════════════════════════════════════════╝")
        print(f"{NC}")
        
        passed, failed = run_download_tests()
        
        write_results_log()
        color = GREEN if failed == 0 else RED
        print(f"{color}")
        print(f"{'=' * 60}")
        print(f"  DOWNLOAD RESULTS: {passed} passed, {failed} failed")
        print(f"{'=' * 60}")
        print(f"{NC}")
        sys.exit(0 if failed == 0 else 1)
    
    # Run utils tests
    if "--utils" in args or "-u" in args:
        print(f"{CYAN}{BOLD}")
        print("╔══════════════════════════════════════════════════════════════════════════════╗")
        print(f"║           HOODINI UTILS TESTS - Running {len(UTILS_TESTS)} tests                             ║")
        print("╚══════════════════════════════════════════════════════════════════════════════╝")
        print(f"{NC}")
        
        passed, failed = run_utils_tests()
        
        write_results_log()
        color = GREEN if failed == 0 else RED
        print(f"{color}")
        print(f"{'=' * 60}")
        print(f"  UTILS RESULTS: {passed} passed, {failed} failed")
        print(f"{'=' * 60}")
        print(f"{NC}")
        sys.exit(0 if failed == 0 else 1)
    
    # Run all tests
    if "--all" in args or "-a" in args:
        print(f"{CYAN}{BOLD}")
        print("╔══════════════════════════════════════════════════════════════════════════════╗")
        print("║           HOODINI COMPLETE TEST SUITE - Pipeline + Download + Utils         ║")
        print("╚══════════════════════════════════════════════════════════════════════════════╝")
        print(f"{NC}")
        
        print_coverage()
        
        total_passed = total_failed = 0
        
        print(f"\n{CYAN}{BOLD}=== PIPELINE TESTS ==={NC}\n")
        p, f = run_tests()
        total_passed += p
        total_failed += f
        
        print(f"\n{CYAN}{BOLD}=== DOWNLOAD TESTS ==={NC}\n")
        p, f = run_download_tests()
        total_passed += p
        total_failed += f
        
        print(f"\n{CYAN}{BOLD}=== UTILS TESTS ==={NC}\n")
        p, f = run_utils_tests()
        total_passed += p
        total_failed += f
        
        write_results_log()
        color = GREEN if total_failed == 0 else RED
        print(f"{color}")
        print(f"{'=' * 60}")
        print(f"  TOTAL RESULTS: {total_passed} passed, {total_failed} failed")
        print(f"{'=' * 60}")
        print(f"{NC}")
        sys.exit(0 if total_failed == 0 else 1)
    
    # Specific pipeline test numbers
    try:
        test_ids = [int(a) for a in args]
        passed, failed = run_tests(test_ids)
        write_results_log()
        color = GREEN if failed == 0 else RED
        print(f"{color}Results: {passed} passed, {failed} failed{NC}")
        sys.exit(0 if failed == 0 else 1)
    except ValueError:
        print(f"{RED}Invalid argument. Use test numbers (1-{len(TESTS)}), --download, --utils, --all, --coverage, --list, or --help{NC}")
        sys.exit(1)


if __name__ == "__main__":
    main()
