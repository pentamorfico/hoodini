"""Pipeline orchestration for the hoodini CLI.

Keeps the CLI thin by encapsulating the main workflow in a single callable
that accepts a typed RuntimeConfig.
"""

from __future__ import annotations

import logging

import polars as pl

from hoodini.config import RuntimeConfig
from hoodini.utils.logging_utils import stage_done, stage_header
from hoodini.utils.memory_utils import reset_tracker

log = logging.getLogger(__name__)


def run_pipeline(config: RuntimeConfig) -> None:
    """Execute the hoodini workflow using the provided config.

    Pipeline Stages and File I/O:
    ==============================

    1. INITIALIZATION (initialize_inputs)
        Expects: config.input_path OR config.inputsheet
        Generates:
        - {output}/ (creates output directory)
        Returns: records DataFrame

    2. IPG PARSING (run_ipg)
        Expects: records DataFrame
        Generates: (no files, enriches records with IPG data)
        Returns: enriched records DataFrame

    3. ASSEMBLY PARSING (run_assembly_parser)
        Expects: records DataFrame, assembly_summary.parquet (packaged data)
        Generates:
        - {output}/assembly_list.txt
        - {output}/assembly_folder/{GCA_*}/*.fna, *.gff
        - {output}/all_neigh.tsv
        - {output}/neighborhood/neighborhoods.fasta
        - {output}/temp.gff
        Returns: all_gff, all_prots, all_neigh DataFrames + valid_uids

    4. PROTEIN COMPARISONS (run_protein_links) [if aai_tree or prot_links]
        Expects: all_prots DataFrame
        Generates:
        - {output}/pairwise_aa.tsv
        Returns: pairwise_aa DataFrame

    5. NUCLEOTIDE COMPARISONS (run_pairwise_nt) [if ani_tree or nt_links]
        Expects: all_neigh, all_gff DataFrames
        Generates:
        - {output}/ani_matrix.tsv (if ani_mode)
        - {output}/nt_links.tsv (if nt_links)
        Returns: pairwise_ani, nt_links DataFrames

    6. PROTEIN CLUSTERING (cluster_proteins)
        Expects: all_prots DataFrame
        Generates:
        - {output}/target_prots.fasta
        - {output}/target_prots.aln (if clust_method != 'none')
        Returns: all_prots with fam_cluster column

    7. PROTEOME SIMILARITY (run_proteome_similarity) [if aai_tree]
        Expects: all_prots, pairwise_aa, all_neigh, all_gff
        Generates:
        - {output}/aai_matrix.tsv
        Returns: pairwise_aai DataFrame

    8. TAXONOMY & TREE (parse_taxonomy_and_build_tree)
        Expects: records, all_gff, all_neigh, all_prots DataFrames
        Generates:
        - {output}/tree.nwk
        - {output}/records.csv
        Returns: tree_str, den_data DataFrame

    9. EXTRA ANNOTATIONS (optional tools)
        - domains (run_domain): {output}/domains.tsv
        - blast (run_blast): enriches all_gff
        - padloc (run_padloc): enriches all_prots
        - emapper (run_emapper): enriches all_prots
        - deffinder (run_defensefinder): enriches all_prots
        - cctyper (run_cctyper): {output}/cctyper/, enriches all_prots + all_gff
        - ncrna (run_ncrna): {output}/ncrna/results.txt, results.sto, enriches all_gff
        - genomad (run_genomad): {output}/genomad/, enriches all_gff

    10. VIZ OUTPUTS (write_viz_outputs)
        Expects: all_gff, all_neigh, all_prots, den_data, tree_str, optional extras
        Generates:
            - {output}/hoodini-viz/parquet/*.parquet (gff, hoods, protein_metadata,
            tree_metadata, nucleotide_links, protein_links, domains, domains_metadata,
            ncrna_metadata)
            - {output}/hoodini-viz/tsv/*.txt (corresponding TSV files)
            - {output}/hoodini-viz/tree.nwk
            - {output}/hoodini-viz/hoodini-viz.html (standalone viewer with embedded data)
    """

    # Initialize memory tracking (only enabled in debug mode)
    tracker = reset_tracker(enabled=config.debug)
    tracker.start_monitoring()

    try:
        _run_pipeline_stages(config, tracker)
    finally:
        tracker.stop_monitoring()
        if config.debug:
            tracker.print_summary()


def _run_pipeline_stages(config: RuntimeConfig, tracker) -> None:
    """Internal function containing all pipeline stages."""

    stage_header("Initializing Hoodini", "🚀")

    from hoodini.pipeline.initialize import initialize_inputs

    with tracker.track_stage("Initialization"):
        records = initialize_inputs(
            input_path=config.input_path,
            inputsheet=config.inputsheet,
            output=config.output,
            force=config.force,
            remote_evalue=config.remote_evalue or 1e-5,
            remote_max_targets=config.remote_max_targets or 100,
        )

    stage_done("Initialization complete")

    stage_header("Parsing IPG data", "🔍")
    from hoodini.pipeline.parse_ipg import run_ipg

    with tracker.track_stage("IPG Parsing"):
        records = run_ipg(
            records_df=records,
            cand_mode=config.cand_mode,
        )

    stage_done("IPG parsing complete")

    stage_header("Downloading and parsing assemblies", "📥")
    from hoodini.pipeline.parse_assemblies import run_assembly_parser

    with tracker.track_stage("Assembly Parsing"):
        result = run_assembly_parser(
            records_df=records,
            output_dir=config.output,
            assembly_folder=config.assembly_folder,
            ncrna=config.ncrna,
            cctyper=config.cctyper,
            genomad=config.genomad,
            blast=config.blast,
            apikey=config.apikey,
            max_concurrent_downloads=config.max_concurrent_downloads,
            num_threads=config.num_threads,
            mod=config.mod,
            wn=config.wn,
            sorfs=config.sorfs,
            minwin=config.minwin,
            minwin_type=config.minwin_type,
        )

    records = result["records"]
    all_gff = result["all_gff"]
    all_prots = result["all_prots"]
    all_neigh = result["all_neigh"]
    valid_uids = result["valid_uids"]

    stage_done("Assembly parsing and neighborhood extraction complete")

    # Abort early if nothing was extracted to avoid downstream errors
    if (all_prots.is_empty() if hasattr(all_prots, "is_empty") else True) or (
        all_neigh.is_empty() if hasattr(all_neigh, "is_empty") else True
    ):
        from hoodini.utils.logging_utils import error

        error("No neighborhoods/proteins extracted; stopping before taxonomy/trees.")
        return

    if config.tree_mode == "aai_tree" or config.prot_links:
        stage_header("Running all-vs-all protein comparisons", "🦠")
        from hoodini.pipeline.protein_links import run_protein_links

        with tracker.track_stage("Protein Comparisons"):
            pairwise_aa = run_protein_links(
                output_dir=config.output,
                all_prots=all_prots,
                threads=config.num_threads,
                evalue=1e-5,
            )

        stage_done("All-vs-all protein comparisons complete")
    else:
        pairwise_aa = None

    if config.tree_mode == "ani_tree" or config.nt_links:
        stage_header("Running pairwise nucleotide comparisons", "🦠")
        from hoodini.pipeline.pairwise_nt import run_pairwise_nt

        with tracker.track_stage("Nucleotide Comparisons"):
            pairwise_ani, nt_links = run_pairwise_nt(
                all_neigh=all_neigh,
                all_gff=all_gff,
                output_dir=config.output,
                nt_aln_mode=config.nt_aln_mode,
                ani_mode=config.ani_mode,
                nt_links=bool(config.nt_links),
                threads=config.num_threads,
            )

        stage_done("Pairwise nucleotide comparisons complete")
    else:
        pairwise_ani = None
        nt_links = None

    stage_header("Clustering neighbor proteins", "✨")
    from hoodini.pipeline.cluster_proteins import cluster_proteins

    with tracker.track_stage("Protein Clustering"):
        all_prots = cluster_proteins(
            all_prots,
            output_dir=config.output,
            clust_method=config.clust_method,
            sorfs=config.sorfs,
        )

    if config.sorfs:
        discarded_sorfs = all_prots.filter(
            pl.col("id").str.contains("sORF") & pl.col("fam_cluster").is_null()
        )
        discarded_sorfs = discarded_sorfs.with_columns(("ID=" + pl.col("id")).alias("gff_id"))
        all_prots = all_prots.filter(~pl.col("id").is_in(discarded_sorfs["id"].unique()))
        all_gff = all_gff.filter(~pl.col("attributes").is_in(discarded_sorfs["gff_id"].unique()))

    stage_done("Clustering complete")

    if config.tree_mode == "aai_tree":
        from hoodini.pipeline.proteome_similarity import run_proteome_similarity

        stage_header("Computing proteome similarity", "🔗")
        with tracker.track_stage("Proteome Similarity"):
            pairwise_aai = run_proteome_similarity(
                all_prots=all_prots,
                pairwise_aa=pairwise_aa,
                all_neigh=all_neigh,
                all_gff=all_gff,
                outdir=config.output,
                mode=config.aai_mode,
                pident_min=config.min_pident,
                subset_mode="target_region",
                win=config.wn,
                win_mode=(
                    config.mod if hasattr(config, "mod") and config.mod is not None else "win_nts"
                ),
                num_threads=config.num_threads,
            )

        stage_done("Proteome similarity complete")
    else:
        pairwise_aai = None

    stage_header("Extracting taxonomic information", "🦠")
    from hoodini.pipeline.taxonomy import parse_taxonomy_and_build_tree

    with tracker.track_stage("Taxonomy & Tree"):
        tree_str, den_data = parse_taxonomy_and_build_tree(
            records=records,
            all_gff=all_gff,
            all_neigh=all_neigh,
            all_prots=all_prots,
            output_dir=config.output,
            tree_mode=config.tree_mode,
            tree_file=config.tree_file,
            num_threads=config.num_threads,
            valid_uids=valid_uids,
            aai_mode=config.aai_mode,
            ani_mode=config.ani_mode,
            aai_subset_mode=config.aai_subset_mode,
            nj_algorithm=config.nj_algorithm,
            pairwise_ani=pairwise_ani,
            pairwise_aai=pairwise_aai,
        )

    domains_data = None
    ncrna_data = None

    with tracker.track_stage("Extra Annotations"):
        if config.domains:
            from hoodini.extra_tools.domain import run_domain

            domains_data = run_domain(all_prots, config.output, config.domains, config.num_threads)

        if config.blast:
            from hoodini.extra_tools.blast import run_blast

            blast_data = run_blast(
                all_neigh, config.output, config.blast, config.num_threads, valid_uids
            )
            if blast_data.height > 0:
                gff_df = pl.DataFrame(
                    {
                        "seqid": blast_data["seqid"],
                        "source": "hoodini",
                        "type": "region",
                        "start": blast_data["start"],
                        "end": blast_data["end"],
                        "score": ".",
                        "strand": "+",
                        "phase": ".",
                        "attributes": "ID=" + blast_data["nc_feature"] + ";",
                    }
                )
                all_gff = pl.concat([all_gff, gff_df], how="vertical")

        if config.padloc:
            from hoodini.extra_tools.padloc import run_padloc

            padloc_df = run_padloc(all_gff, all_prots, config.output, config.num_threads)
            if padloc_df.height > 0:
                all_prots = all_prots.join(padloc_df, on="id", how="left")

        if config.emapper:
            from hoodini.extra_tools.emapper import run_emapper

            emapper_df = run_emapper(all_prots, config.output, config.num_threads)

            if emapper_df.height > 0:
                if "description" in emapper_df.columns and "product" in all_prots.columns:
                    # Create a mapping from id to description
                    desc_map = dict(
                        zip(
                            emapper_df["id"].to_list(),
                            emapper_df["description"].to_list(),
                        )
                    )
                    # Fill empty/null products with emapper descriptions
                    all_prots = all_prots.with_columns(
                        pl.when(
                            pl.col("product").is_null()
                            | (pl.col("product").cast(pl.Utf8).str.strip_chars() == "")
                        )
                        .then(pl.col("id").replace_strict(desc_map, default=pl.col("product")))
                        .otherwise(pl.col("product"))
                        .alias("product")
                    )

                if "id" in emapper_df.columns and "id" in all_prots.columns:
                    all_prots = all_prots.join(emapper_df, on="id", how="left")

        if config.deffinder:
            from hoodini.extra_tools.defensefinder import run_defensefinder

            deffinder_df = run_defensefinder(all_gff, all_prots, config.output)
            if deffinder_df.height > 0:
                all_prots = all_prots.join(deffinder_df, on="id", how="left")

        if config.cctyper:
            from hoodini.extra_tools.cctyper import run_cctyper

            cctyper_df, crispr_df = run_cctyper(
                all_gff, all_prots, all_neigh, config.output, config.num_threads, valid_uids
            )
            if cctyper_df.height > 0:
                all_prots = all_prots.join(cctyper_df, on="id", how="left")
            if crispr_df.height > 0:
                gff_df = crispr_df.select(
                    [
                        pl.col("seqid"),
                        pl.lit("hoodini").alias("source"),
                        pl.lit("region").alias("type"),
                        pl.col("start"),
                        pl.col("end"),
                        pl.lit(".").alias("score"),
                        pl.lit(".").alias("strand"),
                        pl.lit(".").alias("phase"),
                        (pl.lit("ID=") + pl.col("nc_feature") + pl.lit(";")).alias("attributes"),
                    ]
                )
                all_gff = pl.concat([all_gff, gff_df], how="vertical")

        if config.ncrna:
            from hoodini.extra_tools.ncrna import run_ncrna

            ncrna_data = run_ncrna(
                all_neigh, den_data, config.output, config.num_threads, valid_uids, config.ncrna
            )
            if ncrna_data.height > 0:
                gff_df = ncrna_data.select(
                    [
                        pl.col("nucid").alias("seqid"),
                        pl.lit("hoodini").alias("source"),
                        pl.lit("ncRNA").alias("type"),
                        pl.min_horizontal([pl.col("start"), pl.col("end")]).alias("start"),
                        pl.max_horizontal([pl.col("start"), pl.col("end")]).alias("end"),
                        pl.lit(".").alias("score"),
                        pl.col("strand_ncrna").alias("strand"),
                        pl.lit(".").alias("phase"),
                        (pl.lit("ID=") + pl.col("nc_feature") + pl.lit(";")).alias("attributes"),
                    ]
                )
                all_gff = pl.concat([all_gff, gff_df], how="vertical")

        if config.genomad:
            from hoodini.extra_tools.genomad import run_genomad

            genomad_df = run_genomad(all_neigh, config.output, config.num_threads, valid_uids)
            if genomad_df.height > 0:
                gff_df = genomad_df.select(
                    [
                        pl.col("seqid"),
                        pl.lit("hoodini").alias("source"),
                        pl.lit("region").alias("type"),
                        pl.min_horizontal([pl.col("start"), pl.col("end")]).alias("start"),
                        pl.max_horizontal([pl.col("start"), pl.col("end")]).alias("end"),
                        pl.lit(".").alias("score"),
                        pl.lit(".Z").alias("strand"),
                        pl.lit(".").alias("phase"),
                        (pl.lit("ID=") + pl.col("mge_type") + pl.lit(";")).alias("attributes"),
                    ]
                )
                all_gff = pl.concat([all_gff, gff_df], how="vertical")

    stage_done("Extra annotation complete")

    from hoodini.pipeline.write_data import write_viz_outputs

    with tracker.track_stage("Write Outputs"):
        write_viz_outputs(
            output_dir=config.output,
            all_gff=all_gff,
            all_neigh=all_neigh,
            all_prots=all_prots,
            den_data=den_data,
            tree_str=tree_str,
            nt_links=nt_links,
            pairwise_aa=pairwise_aa,
            domains_data=domains_data,
            write_domains=bool(config.domains),
            ncrna_data=ncrna_data,
        )
