import contextlib
from pathlib import Path

import gb_io
import orfipy_core
import polars as pl
import pyrodigal
from Bio.Seq import Seq

from hoodini.utils.logging_utils import info
from hoodini.utils.seq_io import read_fasta


def calculate_overlap(start1, end1, start2, end2):
    """Calculate overlap percentage between two intervals."""
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_start < overlap_end:
        overlap_length = overlap_end - overlap_start
        total_length = max(end1 - start1, end2 - start2)
        return (overlap_length / total_length) * 100
    return 0


def process_features(features, record_version):
    """Extract feature data from GenBank features."""
    feature_data = []
    for feature in features:
        if feature.type == "CDS":
            feature_dict = {
                "seqid": record_version,
                "source": "GenBank",
                "type": "CDS",
                "start": int(feature.location.start),
                "end": int(feature.location.end),
                "strand": "+" if feature.location.strand == 1 else "-",
                "phase": ".",
            }
            # Extract qualifiers
            for qual in ["protein_id", "product", "gene"]:
                if qual in feature.qualifiers:
                    feature_dict[qual] = feature.qualifiers[qual][0]
            feature_data.append(feature_dict)
    return feature_data


def unwrap_attributes(gff_df):
    """Unwrap attributes column in GFF DataFrame."""
    if "attributes" in gff_df.columns:
        # Parse attributes into separate columns
        attributes_list = []
        for attr in gff_df["attributes"]:
            attr_dict = {}
            if attr and isinstance(attr, str):
                for pair in attr.split(";"):
                    if "=" in pair:
                        key, value = pair.split("=", 1)
                        attr_dict[key.strip()] = value.strip()
            attributes_list.append(attr_dict)
        attributes_df = pl.DataFrame(attributes_list)
        gff_df = pl.concat([gff_df.drop("attributes"), attributes_df], how="horizontal")
    return gff_df


def extract_neighborhood(
    protein_id,
    nucleotide_id,
    gbf_file,
    gff_file,
    faa_file,
    fna_file,
    mode="win_nts",
    window=None,
    strand=None,
    start=None,
    end=None,
    unique_id=None,
    input_type=None,
    sorfs=None,
):
    info(f"✔️\tExtracting neighborhood {unique_id}")
    neighborhood = {}
    if gbf_file:

        record_found = False

        if not Path(gbf_file).exists():
            return (
                None,
                None,
                unique_id,
                "GenBank file not found",
            )

        try:
            record_iter = gb_io.iter(gbf_file)
        except Exception:
            return None, None, unique_id, "gb_io failed to open GenBank iterator"

        if nucleotide_id:
            record_found = False
            for record in record_iter:
                record_version = getattr(record, "version", None)
                if record_version and nucleotide_id in record_version:
                    record_found = True
                    break
        if record_found:
            feature_data = process_features(record.features, record_version)
            feature_data = pl.DataFrame(feature_data)
            if "attributes" in feature_data.columns:
                attributes_df = pl.DataFrame(feature_data["attributes"].to_list())
                attributes_df.drop(["protein_id"])
                feature_data = pl.concat(
                    [feature_data.drop(["attributes"]), attributes_df], how="horizontal"
                )
                feature_data = feature_data.rename({"translation": "sequence"})
            else:
                return None, None, unique_id, "GenBank file is not annotated"
        else:
            return None, None, unique_id, "GenBank record not found"

        if input_type == "protein":

            if "protein_id" in feature_data.columns and not (start and end):
                start = feature_data[feature_data["protein_id"] == protein_id]["start"].iloc[0]
                end = feature_data[feature_data["protein_id"] == protein_id]["end"].iloc[0]
                strand = feature_data[feature_data["protein_id"] == protein_id]["strand"].iloc[0]
            start, end = int(start), int(end)

            if mode == "win_nts":
                start_win = start - window
                end_win = end + window
                start_win = max(start_win, 0)
                end_win = min(end_win, len(record.sequence))
                subgff = feature_data.query("start>=@start_win & end<=@end_win")

            elif mode == "win_ngen":
                subgff = feature_data.reset_index(drop=True)
                prot_index = subgff[subgff["protein_id"] == protein_id].index.to_list()[0]
                subgff = subgff[prot_index - window : prot_index + window]
                start_win = subgff["start"].min()
                end_win = subgff["end"].max()

        elif input_type == "nucleotide":

            if nucleotide_id and (start and end) and window:
                start = int(start)
                end = int(end)
                if mode == "win_nts":
                    start_win = start - window
                    end_win = end + window
                    if not strand:
                        strand = "+" if end > start else "-"
                    subgff = feature_data.query("start>=@start_win & end<=@end_win")

                elif mode == "win_ngen":
                    subgff = feature_data.reset_index(drop=True)
                    start_index = subgff[subgff["start"] >= start].index.to_list()[0]
                    end_index = subgff[subgff["end"] <= end].index.to_list()[-1]
                    subgff = subgff[start_index - window : end_index + window]
                    start_win = subgff["start"].min()
                    end_win = subgff["end"].max()
                    if not strand:
                        strand = "+" if end > start else "-"

            elif not window and nucleotide_id and (start and end):
                start, end = int(start), int(end)
                if not strand:
                    strand = "-" if end < start else "+"
                subgff = feature_data.query(
                    "seqid == @nucleotide_id & type =='CDS' & start>=@start & end<=@end"
                )
                start_win = subgff["start"].min()
                end_win = subgff["end"].max()

            elif not (start and end):
                start = end = 0
                subgff = feature_data
                start_win = subgff["start"].min()
                end_win = subgff["end"].max()
                strand = "+"

        start_win = max(start_win, 0)
        end_win = min(end_win, len(record.sequence))
        subgff["id"] = subgff["protein_id"]
        header = [
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "protein_id",
            "id",
            "sequence",
        ]
        if "product" in subgff.columns:
            header.append("product")
        else:
            subgff["product"] = None
        subgff = subgff[header]
        neighborhood = {
            "seqid": record_version,
            "start_target": start,
            "end_target": end,
            "start_win": start_win,
            "end_win": end_win,
            "strand_win": strand,
            "sequence": record.sequence[start_win:end_win].decode("utf-8"),
            "unique_id": unique_id,
        }
        if sorfs:
            orf_finder = pyrodigal.GeneFinder(meta=True, min_gene=10, max_overlap=9)
            new_genes = []

            for i, pred in enumerate(
                orf_finder.find_genes(record.sequence[start_win:end_win].decode("utf-8"))
            ):
                overlap_flag = False
                for row in subgff.iter_rows(named=True):
                    overlap_percentage = calculate_overlap(
                        row["start"], row["end"], pred.begin + start_win, pred.end + start_win
                    )
                    if overlap_percentage > 10:
                        overlap_flag = True
                        break
                if not overlap_flag:
                    new_genes.append(
                        {
                            "seqid": nucleotide_id,
                            "source": "pyrodigal",
                            "type": "CDS",
                            "start": pred.begin + start_win,
                            "end": pred.end + start_win,
                            "score": pred.score,
                            "strand": "-" if pred.strand == "-1" else "+",
                            "phase": ".",
                            "protein_id": f"sORF_{unique_id}_{i}",
                            "id": f"sORF_{unique_id}_{i}",
                            "sequence": pred.translate(),
                        }
                    )

            if new_genes:
                new_genes_df = pl.DataFrame(new_genes)
                subgff = pl.concat([subgff, new_genes_df], how="vertical")

            new_genes = []
            seq = record.sequence[start_win:end_win].decode("utf-8").upper()
            for i, (start, stop, strand, _description) in enumerate(
                orfipy_core.orfs(seq, minlen=100, maxlen=1000, partial3=False, between_stops=False)
            ):
                overlap_flag = False
                for row in subgff.iter_rows(named=True):
                    overlap_percentage = calculate_overlap(
                        row["start"], row["end"], start + start_win, stop + start_win
                    )
                    if overlap_percentage > 0:
                        overlap_flag = True
                        break

                if not overlap_flag:
                    orf_sequence = Seq(record.sequence[start_win:end_win][start:stop])
                    if strand == "-":
                        orf_sequence = orf_sequence.reverse_complement()
                    protein_sequence = orf_sequence.translate(table=11, to_stop=True)

                    new_genes.append(
                        {
                            "seqid": nucleotide_id,
                            "source": "orfipy",
                            "type": "CDS",
                            "start": start + start_win,
                            "end": stop + start_win,
                            "score": ".",
                            "strand": "-" if strand == "-" else "+",
                            "phase": ".",
                            "protein_id": f"sORF_orfipy_{unique_id}_{i}",
                            "id": f"sORF_orfipy_{unique_id}_{i}",
                            "sequence": protein_sequence,
                        }
                    )

            if new_genes:
                new_genes_df = pl.DataFrame(new_genes)
                subgff = pl.concat([subgff, new_genes_df], how="vertical")

        neighborhood = pl.DataFrame(neighborhood, index=[0])

    elif gff_file and faa_file:
        gff_header = [
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ]
        if not Path(gff_file).exists():
            return (
                None,
                None,
                unique_id,
                "GFF file not found",
            )
        if not Path(faa_file).exists():
            return None, None, unique_id, "FAA file not found"

        try:
            gff = pl.read_csv(
                filepath_or_buffer=gff_file,
                separator="\t",
                comment="#",
                names=gff_header,
                engine="c",
            )
        except Exception:
            return None, None, unique_id, "Failed to read GFF file"
        try:
            faa_df = read_fasta(faa_file)
        except Exception:
            return None, None, unique_id, "Failed to read FAA file"

        if input_type == "protein":
            if protein_id and window:
                query = f"={protein_id}"
                start = gff[gff["attributes"].str.contains(query)]["start"].to_list()
                if not start:
                    return None, None, unique_id, "Protein not found in GFF file"
                else:
                    start = start[0]
                end = gff[gff["attributes"].str.contains(query)]["end"].to_list()[0]
                strand = gff[gff["attributes"].str.contains(query)]["strand"].to_list()[0]
                nucleotide_id = gff[gff["attributes"].str.contains(query)]["seqid"].to_list()[0]

                if mode == "win_nts":
                    start_win = start - window
                    end_win = end + window
                    start_win = max(start_win, 0)
                    gff_nuc = gff.query("seqid == @nucleotide_id")
                    end_win = min(end_win, gff_nuc["end"].max())
                    subgff = gff.query(
                        "seqid == @nucleotide_id & type =='CDS' & start>=@start_win & end<=@end_win"
                    )
                elif mode == "win_ngen":
                    subgff = gff.query("seqid == @nucleotide_id & type =='CDS'").reset_index(
                        drop=True
                    )
                    prot_index = subgff[subgff["attributes"].str.contains(query)].index.to_list()[0]
                    subgff = subgff[prot_index - window : prot_index + window]
                    start_win = subgff["start"].min()
                    end_win = subgff["end"].max()

                if strand == "-":
                    pass

        elif input_type == "nucleotide":
            if nucleotide_id and (start and end) and window:
                start, end = int(start), int(end)
                if not strand:
                    strand = "-" if end < start else "+"

                if mode == "win_nts":
                    start_win = start - window
                    end_win = end + window
                    start_win = max(start_win, 0)
                    subgff = gff.query(
                        "seqid == @nucleotide_id & type =='CDS' & start>=@start_win & end<=@end_win"
                    )
                elif mode == "win_ngen":
                    subgff = gff.query("seqid == @nucleotide_id & type =='CDS'").reset_index(
                        drop=True
                    )
                    prot_index = subgff[subgff["attributes"].str.contains(query)].index.to_list()[0]
                    subgff = subgff[prot_index - window : prot_index + window]
                    start_win = subgff["start"].min()
                    end_win = subgff["end"].max()

                if strand == "-":
                    pass

            elif not window and nucleotide_id and (start and end):
                start, end = int(start), int(end)
                if not strand:
                    strand = "-" if end < start else "+"
                subgff = gff.query(
                    "seqid == @nucleotide_id & type =='CDS' & start>=@start & end<=@end"
                )
                start_win = start
                end_win = end
                if strand == "-":
                    pass

            elif (
                not window
                and nucleotide_id
                and not (start and end)
                or window
                and nucleotide_id
                and not (start and end)
            ):
                start = end = 0
                subgff = gff.query("seqid == @nucleotide_id & type =='CDS'")
                start_win = 0
                end_win = subgff["end"].max()
                strand = "+"

            else:
                return (
                    None,
                    None,
                    unique_id,
                    "Invalid usage of parameters",
                )

        subgff = unwrap_attributes(subgff)
        key_join = "protein_id" if "protein_id" in subgff.columns else "ID"

        subgff = subgff.join(
            faa_df[["id", "sequence"]], left_on=key_join, right_on="id", how="left"
        )

        if fna_file:
            fna_df = read_fasta(fna_file)
            nucleotide_id = str(nucleotide_id)
            faa_df["id"] = faa_df["id"].astype(str)
            if nucleotide_id in fna_df["id"].to_list():
                sequence = fna_df[fna_df["id"] == nucleotide_id]["sequence"].to_list()[0]
                end_win = end + window
                end_win = min(end_win, len(sequence))
                if sorfs:
                    orf_finder = pyrodigal.GeneFinder(meta=True)
                    new_genes = []

                    for i, pred in enumerate(orf_finder.find_genes(sequence.encode())):
                        overlap_flag = False
                        for row in subgff.iter_rows(named=True):
                            overlap_percentage = calculate_overlap(
                                row["start"], row["end"], pred.begin, pred.end
                            )
                            if overlap_percentage > 5:
                                overlap_flag = True
                                break

                        if not overlap_flag:
                            new_genes.append(
                                {
                                    "seqid": nucleotide_id,
                                    "source": "pyrodigal",
                                    "type": "CDS",
                                    "start": pred.begin + start_win,
                                    "end": pred.end + start_win,
                                    "score": pred.score,
                                    "strand": "-" if pred.strand == "-1" else "+",
                                    "phase": ".",
                                    key_join: f"{key_join}=sORF_{unique_id}_{i}",
                                    "sequence": pred.translate(),
                                }
                            )

                    if new_genes:
                        new_genes_df = pl.DataFrame(new_genes)
                        subgff = pl.concat([subgff, new_genes_df], how="vertical")
        else:
            sequence = None
        neighborhood = {
            "seqid": nucleotide_id,
            "start_target": start,
            "end_target": end,
            "start_win": start_win,
            "end_win": end_win,
            "strand_win": strand,
            "sequence": sequence[start_win:end_win],
            "unique_id": unique_id,
        }
        neighborhood = pl.DataFrame(neighborhood, index=[0])

    subgff["target_prot"] = protein_id
    subgff["target_nuc"] = nucleotide_id
    subgff["unique_id"] = str(unique_id)

    try:
        if isinstance(subgff, pl.DataFrame):
            if "id" not in subgff.columns:
                for cand in ("protein_id", "ID", "gene_id"):
                    if cand in subgff.columns:
                        subgff["id"] = subgff[cand]
                        break

            for redundant in ("ID", "gene_id"):
                if redundant in subgff.columns and redundant != "id":
                    with contextlib.suppress(Exception):
                        subgff.drop([redundant])

            if "id" in subgff.columns:
                with contextlib.suppress(Exception):
                    subgff["id"] = subgff["id"].astype(str)
    except Exception:
        pass

    if "product" not in subgff.columns:
        subgff["product"] = None

    return subgff, neighborhood, unique_id
