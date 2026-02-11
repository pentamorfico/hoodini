"""Runtime configuration helpers.

Centralizes how defaults, user TOML files, and CLI overrides are merged into a
single typed object. Keep this layer free of CLI concerns so it can be reused
by tests or other entry points.
"""

from __future__ import annotations

from collections.abc import Mapping, MutableMapping
from dataclasses import dataclass, field, replace
from typing import Any


def _lower_keys(data: Mapping[str, Any]) -> dict[str, Any]:
    """Return a new dict with string keys lowercased (TOML keys often vary)."""
    lowered: dict[str, Any] = {}
    for key, value in data.items():
        lowered[key.lower()] = value
    return lowered


def _flatten_config(grouped: Mapping[str, Any]) -> dict[str, Any]:
    """Flatten a nested TOML-style mapping into a single-level dict."""
    flat: dict[str, Any] = {}
    for section, values in grouped.items():
        if isinstance(values, MutableMapping):
            flat.update(_lower_keys(values))
        else:
            flat[section] = values
    return flat


@dataclass(slots=True)
class RuntimeConfig:
    input_path: str | None = None
    inputsheet: str | None = None
    output: str | None = None

    max_concurrent_downloads: int | None = None
    apikey: str | None = None
    num_threads: int | None = None
    assembly_folder: str | None = None

    prot_links: bool = False
    nt_links: bool = False

    ani_mode: str | None = None
    nt_aln_mode: str | None = None
    blast: str | None = None
    cand_mode: str | None = None
    clust_method: str | None = None
    mod: str | None = None
    wn: int | None = None
    minwin: int | None = None
    minwin_type: str | None = None

    tree_mode: str | None = None
    tree_file: str | None = None
    aai_mode: str | None = None
    aai_subset_mode: str | None = None
    nj_algorithm: str | None = None
    remote_evalue: float | None = None
    remote_max_targets: int | None = None

    padloc: bool = False
    deffinder: bool = False
    ncrna: str | None = None
    cctyper: bool = False
    trna: bool = False
    genomad: bool = False
    sorfs: bool = False
    domains: list[str] | None = field(default=None)
    emapper: bool = False
    min_pident: float = 30.0

    keep: bool = False
    force: bool = False
    debug: bool = False

    def replace(self, **kwargs: Any) -> RuntimeConfig:
        """Return a copy with provided fields updated."""
        return replace(self, **kwargs)


def build_runtime_config(
    *,
    defaults: Mapping[str, Any],
    file_overrides: Mapping[str, Any] | None = None,
    cli_overrides: Mapping[str, Any] | None = None,
) -> RuntimeConfig:
    """Merge defaults + file + CLI into a RuntimeConfig.

    Later sources override earlier ones. Unknown keys are ignored to keep the
    dataclass strict.
    """

    merged: dict[str, Any] = {}
    for source in (defaults, file_overrides or {}, cli_overrides or {}):
        merged.update(_flatten_config(source))

    if "aai-subset-mode" in merged and "aai_subset_mode" not in merged:
        merged["aai_subset_mode"] = merged["aai-subset-mode"]
    if "aai_subset_mode" in merged and "aai-subset-mode" not in merged:
        merged["aai-subset-mode"] = merged["aai_subset_mode"]

    allowed_fields: set[str] = {field.name for field in RuntimeConfig.__dataclass_fields__.values()}  # type: ignore[attr-defined]
    filtered: dict[str, Any] = {k: v for k, v in merged.items() if k in allowed_fields}

    config = RuntimeConfig(**filtered)
    if config.input_path is None:
        config.input_path = None
    if config.inputsheet is None:
        config.inputsheet = None
    return config


__all__ = ["RuntimeConfig", "build_runtime_config"]
