"""Small collection helpers."""

from __future__ import annotations


def flat(seq):
    return [item for sublist in seq for item in sublist]
