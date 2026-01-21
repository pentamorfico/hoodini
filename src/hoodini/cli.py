import logging
import os
import sys
import warnings
from pathlib import Path

import rich_click as click
import tomli
from click.core import ParameterSource
from rich.logging import RichHandler

from hoodini.config import build_runtime_config, load_default_config
from hoodini.pipeline.runner import run_pipeline
from hoodini.utils.cli_helpers import MutuallyExclusiveOption
from hoodini.utils.logging_utils import (
    configure_logging,
    console,
    header,
    stage_done,
    stage_header,
    success,
)
from hoodini.utils.validation import validate_domains, validate_input_file

click.rich_click.USE_RICH_MARKUP = True
click.rich_click.STYLE_ERRORS_SUGGESTIONS = "yellow"
click.rich_click.STYLE_HELP_OPTIONS = "bold cyan"
click.rich_click.STYLE_HELP_OPTIONS_DEFAULTS = "dim"

# Group CLI options by theme
click.rich_click.OPTION_GROUPS = {
    "hoodini run": [
        {
            "name": "Input/Output",
            "options": ["--config", "--input", "--inputsheet", "--output", "--force", "--keep"],
        },
        {
            "name": "Performance",
            "options": ["--max-concurrent-downloads", "--num-threads", "--api-key"],
        },
        {
            "name": "Data Sources",
            "options": ["--assembly-folder", "--blast"],
        },
        {
            "name": "Neighborhood Window",
            "options": ["--win-mode", "--win", "--min-win", "--min-win-type"],
        },
        {
            "name": "Clustering",
            "options": ["--cand-mode", "--clust-method"],
        },
        {
            "name": "Tree Construction",
            "options": ["--tree-mode", "--tree-file"],
        },
        {
            "name": "Pairwise Comparisons",
            "options": [
                "--prot-links",
                "--nt-links",
                "--ani-mode",
                "--nt-aln-mode",
                "--aai-mode",
                "--aai-subset-mode",
                "--min-pident",
            ],
        },
        {
            "name": "Remote BLAST",
            "options": ["--remote-evalue", "--remote-max-targets"],
        },
        {
            "name": "Annotations",
            "options": [
                "--padloc",
                "--deffinder",
                "--cctyper",
                "--ncrna",
                "--genomad",
                "--sorfs",
                "--emapper",
                "--domains",
            ],
        },
        {
            "name": "Logging",
            "options": ["--quiet", "--debug"],
        },
    ],
}


@click.group()
def cli():
    """🦉 hoodini: gene-centric comparative genomic analysis using publicly available data"""


@cli.command()
@click.option(
    "--config",
    "config_file",
    default=None,
    type=click.Path(exists=True),
    help="TOML config file to load parameters from.",
)
@click.option(
    "--input",
    "input_path",
    cls=MutuallyExclusiveOption,
    mutually_exclusive=["inputsheet"],
    default=None,
    type=str,
    callback=validate_input_file,
    help="Path to a single-column input file or a literal protein ID/FASTA (mutually exclusive with --inputsheet).",
)
@click.option(
    "--inputsheet",
    cls=MutuallyExclusiveOption,
    mutually_exclusive=["input_path"],
    default=None,
    type=click.Path(exists=True),
    callback=validate_input_file,
    help="Path to a TSV input file (mutually exclusive with --input).",
)
@click.option("--output", help="Output folder name.")
@click.option(
    "--max-concurrent-downloads",
    "max_concurrent_downloads",
    type=int,
    help="Maximum concurrent downloads.",
)
@click.option("--api-key", "apikey", help="NCBI API key.")
@click.option("--num-threads", "num_threads", type=int, help="Number of threads.")
@click.option("--assembly-folder", "assembly_folder", help="Path to a local assembly folder.")
@click.option("--prot-links", "prot_links", is_flag=True, help="Run pairwise protein comparisons.")
@click.option("--nt-links", "nt_links", is_flag=True, help="Run pairwise nucleotide comparisons.")
@click.option(
    "--ani-mode",
    "ani_mode",
    type=click.Choice(["skani", "blastn"]),
    help="Choose ANI calculation method for ANI trees.",
)
@click.option(
    "--nt-aln-mode",
    "nt_aln_mode",
    type=click.Choice(["blastn", "fastani", "minimap2", "intergenic_blastn"]),
    help="Nucleotide alignment mode for pairwise comparisons.",
)
@click.option("--blast", help="BLAST query file to use.")
@click.option(
    "--cand-mode",
    "cand_mode",
    type=click.Choice(["any_ipg", "best_ipg", "best_id", "one_id", "same_id"]),
    help="Mode for selecting IPG representative.",
)
@click.option(
    "--clust-method",
    "clust_method",
    type=click.Choice(["diamond_deepclust", "deepmmseqs", "jackhmmer", "blastp"]),
    help="Protein clustering method.",
)
@click.option(
    "--win-mode",
    "mod",
    type=click.Choice(["win_nts", "win_genes"]),
    help="Window mode: 'win_nts' or 'win_genes'.",
)
@click.option("--win", "wn", type=int, help="Window size (genes or nucleotides).")
@click.option("--min-win", "minwin", type=int, help="Min window size on each side.")
@click.option(
    "--min-win-type",
    "minwin_type",
    type=click.Choice(["total", "upstream", "downstream", "both"]),
    help="Type of min window.",
)
@click.option(
    "--tree-mode",
    "tree_mode",
    type=click.Choice(
        [
            "taxonomy",
            "fast_nj",
            "aai_tree",
            "ani_tree",
            "fast_ml",
            "use_input_tree",
            "foldmason_tree",
            "neigh_similarity_tree",
            "neigh_phylo_tree",
        ]
    ),
    help="Tree building method.",
)
@click.option("--tree-file", "tree_file", help="Path to the tree file.")
@click.option(
    "--aai-mode",
    "aai_mode",
    type=click.Choice(["wgrr", "aai"]),
    help="Mode for AAI tree construction (wgrr or aai). Note: 'hyper' mode is not supported for AAI trees.",
)
@click.option(
    "--aai-subset-mode",
    "aai_subset_mode",
    type=click.Choice(["target_prot", "target_region", "window"]),
    help="Subset mode for selecting proteins in AAI tree construction.",
)
@click.option(
    "--remote-evalue",
    "remote_evalue",
    type=float,
    help="E-value for remote BLAST when providing a single protein ID/FASTA as input.",
)
@click.option(
    "--remote-max-targets",
    "remote_max_targets",
    type=int,
    help="Maximum targets to retrieve in remote BLAST for single protein input.",
)
@click.option("--padloc", is_flag=True, help="Run PADLOC for antiphage defense.")
@click.option("--deffinder", is_flag=True, help="Run DefenseFinder for antiphage defense.")
@click.option(
    "--ncrna",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Path to CM models file for ncRNA prediction with Infernal (e.g., Rfam.cm).",
)
@click.option("--cctyper", is_flag=True, help="Run CCtyper for CRISPR-Cas prediction.")
@click.option("--genomad", is_flag=True, help="Run GenoMAD for MGE identification.")
@click.option("--sorfs", is_flag=True, help="Reannotate small open reading frames.")
@click.option(
    "--domains",
    callback=validate_domains,
    help="Comma-separated list of MetaCerberus domain databases.",
)
@click.option(
    "--emapper",
    "emapper",
    is_flag=True,
    default=False,
    help="Run eggNOG-mapper (emapper) to annotate proteins and append annotations to protein metadata.",
)
@click.option(
    "--min-pident",
    "min_pident",
    type=float,
    default=30.0,
    help="Minimum percent identity threshold for BLAST hits in wGRR/AAI calculations.",
)
@click.option("--keep", is_flag=True, help="Keep temporary files (do not delete).")
@click.option("--force", is_flag=True, help="Overwrite existing output folder if it exists.")
@click.option("--quiet", is_flag=True, help="Silence all non-error output.")
@click.option("--debug", is_flag=True, help="Enable verbose debug logging.")
@click.pass_context
def run(ctx, config_file: str | None, quiet: bool, debug: bool, **cli_kwargs) -> None:
    """
    Run hoodini with default parameters or from a config file.
    """
    log_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(level=log_level, format="%(message)s", handlers=[RichHandler()])
    if not debug:
        logging.getLogger("urllib3").setLevel(logging.WARNING)
        logging.getLogger("requests").setLevel(logging.WARNING)
        logging.getLogger("httpx").setLevel(logging.WARNING)
        warnings.filterwarnings(
            "ignore",
            message="pkg_resources is deprecated as an API.*",
            category=UserWarning,
        )
    configure_logging(quiet=quiet, debug=debug)
    subtitle = "large-scale gene neighborhood analyses that feel like magic"
    header("hoodini 🦉🎩", subtitle, border_style="light_slate_grey")
    # Separator just under the top header; width ~ text width
    sep_width = max(len("hoodini 🦉🎩"), len(subtitle)) + 4
    bar = "━" * sep_width
    console.print(f"[light_slate_grey]{bar}[/light_slate_grey]")

    file_kwargs = {}
    if config_file:
        with open(config_file, "rb") as f:
            grouped = tomli.load(f)
            file_kwargs = {k: v for section in grouped.values() for k, v in section.items()}

    cli_clean = {}
    for k, v in cli_kwargs.items():
        src = ctx.get_parameter_source(k)
        if src in (ParameterSource.DEFAULT, ParameterSource.DEFAULT_MAP):
            continue
        cli_clean[k] = v

    config = build_runtime_config(
        defaults=load_default_config(),
        file_overrides=file_kwargs,
        cli_overrides={**cli_clean, "debug": debug},
    )

    config = config.replace(
        ani_mode=config.ani_mode if config.tree_mode == "ani_tree" else None,
        nt_aln_mode=(
            config.nt_aln_mode if (config.nt_links or config.tree_mode == "ani_tree") else None
        ),
    )

    if not config.input_path and not config.inputsheet:
        raise click.UsageError("One of '--input' or '--inputsheet' must be provided.")

    ctx.ensure_object(dict)
    ctx.obj["config"] = config

    run_pipeline(config)


@cli.group()
def download():
    """Download resources used by Hoodini."""


@download.command("assembly_summary")
def download_assembly_summary():
    """Download or update the assembly_summary.parquet database."""
    from hoodini.download.assembly_summary import download_assembly_summary_db

    stage_header("Downloading NCBI Assembly Database", "📥")
    asm_summary_path = download_assembly_summary_db()
    stage_done(f"Saved to {asm_summary_path}")


@download.command("metacerberus")
@click.argument("dbs", required=False, default="all")
@click.option("--force", is_flag=True, help="Overwrite existing files.")
def download_metacerberus(dbs, force):
    """Download MetaCerberus HMM/TSV databases from OSF.io."""
    from hoodini.download.metacerberus import main as metacerberus_main

    if not dbs or dbs.strip().lower() == "all":
        metacerberus_main(None, force=force)
    else:
        metacerberus_main(dbs, force=force)
        stage_done(f"MetaCerberus download of {dbs} complete!")


@download.command("type_dive")
def download_type_dive():
    """Download and normalize BacDive and PhageDive DSMZ databases."""
    from hoodini.download.type_dive import main as type_dive_main

    stage_header("Downloading DSMZ BacDive/PhageDive databases", "🦠")
    type_dive_main()
    stage_done("DSMZ BacDive/PhageDive download and normalization complete!")


@download.command("contig_lengths")
@click.option(
    "--api-key",
    "api_key",
    default=os.environ.get("NCBI_API_KEY"),
    help="NCBI API key (overrides environment variable NCBI_API_KEY).",
)
@click.option(
    "--skip-assembly-summary",
    "skip_assembly_summary",
    is_flag=True,
    default=False,
    help="Skip refreshing local assembly_summary.parquet (use existing local copy).",
)
def download_contig_lengths(api_key, skip_assembly_summary):
    """Download missing NCBI contig length records and update precomputed list."""
    stage_header("Downloading NCBI contig lengths", "📥")
    from hoodini.download.contig_lengths import download_contig_lengths as impl

    impl(api_key=api_key, skip_assembly_summary=skip_assembly_summary)
    stage_done("NCBI contig length download complete")


@download.command("databases")
@click.option("--force", is_flag=True, help="Force re-download of files")
@click.option("--skip-padloc", is_flag=True, help="Skip padloc DB update.")
@click.option("--skip-deffinder", is_flag=True, help="Skip defense-finder model install.")
@click.option("--skip-genomad", is_flag=True, help="Skip GenoMAD download.")
@click.option("--skip-emapper", is_flag=True, help="Skip downloading emapper/mmseqs DB.")
@click.option("--skip-parquet", is_flag=True, help="Skip downloading eggNOG parquet support files.")
@click.option(
    "--skip-contig-lengths", is_flag=True, help="Skip downloading contig_lengths.parquet."
)
@click.option(
    "--skip-typedive", is_flag=True, help="Skip downloading BacDive/PhageDive databases."
)
@click.option(
    "--threads",
    "num_threads",
    type=int,
    default=0,
    help="Number of threads for aria2c and pigz (0 = use all cores).",
)
def download_databases(
    force,
    skip_padloc,
    skip_deffinder,
    skip_genomad,
    skip_emapper,
    skip_parquet,
    skip_contig_lengths,
    skip_typedive,
    num_threads,
):
    """Download support databases used by some extra tools (emapper, PADLOC models, etc.)"""
    from hoodini.download.databases import main as db_main

    db_main(
        force=force,
        skip_padloc=skip_padloc,
        skip_deffinder=skip_deffinder,
        skip_genomad=skip_genomad,
        skip_emapper=skip_emapper,
        skip_parquet=skip_parquet,
        skip_contig_lengths=skip_contig_lengths,
        skip_typedive=skip_typedive,
        num_threads=num_threads,
    )


@cli.group()
def utils():
    """Utility commands for Hoodini (e.g., sequence metadata helpers)."""


@utils.command("nuc2asmlen")
@click.argument("input_file", type=click.Path(exists=True))
@click.option(
    "--output",
    "output_file",
    type=click.Path(),
    default=None,
    help="Optional output file (TSV). If not set, prints to stdout.",
)
def nuc2asmlen(input_file, output_file):
    """Fetch assembly and length metadata for nuccore/contig accessions (with fallback)."""
    from hoodini.pipeline.helpers.nuc2asmlen import run_nuc2asmlen

    df = run_nuc2asmlen(input_file)

    if output_file:
        df.write_csv(output_file, separator="\t")
        success(f"Saved results to {output_file}")
    else:
        sys.stdout.write(df.write_csv(separator="\t", include_header=False))


@utils.command("prefetch_links")
@click.argument("input_file", type=click.Path(exists=True))
@click.option(
    "--output",
    "output_file",
    type=click.Path(),
    default=None,
    help="Optional output file (TSV). If not set, prints to stdout.",
)
@click.option(
    "--kinds",
    "kinds",
    default="gbff,gff,fna,faa,sequence_report",
    help="Comma-separated file kinds to generate links for.",
)
def prefetch_links(input_file, output_file, kinds):
    """Generate prefetched links (assembly_id, file_type, link) for a list of assemblies."""
    from hoodini.pipeline.helpers.prefetch_links import get_prefetched_link_table

    with open(input_file) as fh:
        accs = [l.strip() for l in fh if l.strip()]

    kinds_list = [k.strip() for k in kinds.split(",") if k.strip()]
    df = get_prefetched_link_table(accs, kinds=kinds_list)
    if output_file:
        df.write_csv(
            output_file,
            separator="\t",
            include_header=False,
            columns=["assembly_id", "file_type", "link"],
        )
        success(f"Saved prefetched links to {output_file}")
    else:
        sys.stdout.write(
            df.write_csv(
                separator="\t", include_header=False, columns=["assembly_id", "file_type", "link"]
            )
        )


if __name__ == "__main__":
    cli()

main = cli
