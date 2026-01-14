"""ID parsing and categorization utilities."""

from __future__ import annotations

import re


def is_refseq_nuccore(nuc_id) -> bool:
    """Return True if the nuccore accession is a RefSeq accession else False."""
    refseq_prefixes = ("NC_", "NZ_", "NM_", "NR_", "XM_", "XR_", "AP_", "YP_", "XP_", "WP_")
    return isinstance(nuc_id, str) and nuc_id.startswith(refseq_prefixes)


def switch_assembly_prefix(asm_id):
    if not isinstance(asm_id, str):
        return asm_id
    if asm_id.startswith("GCA_"):
        return "GCF_" + asm_id[4:]
    if asm_id.startswith("GCF_"):
        return "GCA_" + asm_id[4:]
    return asm_id


def categorize_id(id_: str) -> dict[str, str | None]:
    parts = id_.split(":")
    id_part = parts[0]

    uniprot_pattern = re.compile(
        r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|"
        r"[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]"
        r"(?:[A-Z][A-Z0-9]{2}[0-9])?)$"
    )
    nucleotide_patterns = [
        re.compile(
            r"^("
            + "|".join(
                [
                    "NC",
                    "NG",
                    "NM",
                    "NR",
                    "NT",
                    "NW",
                    "NZ",
                    "AC",
                    "AP",
                    "MT",
                    "PP",
                    "OR",
                    "OZ",
                    "LR",
                    "LN",
                    "KX",
                ]
            )
            + r")(_[A-Z]+\d+|\d+)(\.\d+)?(:\d+-\d+)?$"
        ),
        re.compile(r"^[A-Z]{1,2}\d{5,8}(\.\d+)?$"),
        re.compile(r"^[A-Z]{4,6}\d{8,}(\.\d+)?$"),
    ]
    protein_patterns = [
        re.compile(r"^(" + "|".join(["NP", "XP", "YP", "WP", "ZP"]) + r")_\d+(\.\d+)?$"),
        re.compile(r"^[A-Z]{3}\d{5,8}(\.\d+)?$"),
    ]

    if re.match(uniprot_pattern, id_part):
        return {"type": "uniprot", "id": id_part, "protein_id": None}
    if any(re.match(pattern, id_part) for pattern in nucleotide_patterns):
        return {
            "type": "nucleotide",
            "id": id_part,
            "protein_id": parts[1] if len(parts) > 1 else None,
        }
    if any(re.match(pattern, id_part) for pattern in protein_patterns):
        return {"type": "protein", "id": id_part, "protein_id": None}
    return {"type": "unmatched", "id": id_part, "protein_id": None}
