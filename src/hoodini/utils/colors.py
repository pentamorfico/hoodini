"""Color utility helpers."""

from __future__ import annotations


def desaturate(hex_color: str, amount: float = 0.2) -> str:
    hex_color = hex_color.lstrip("#")
    rgb = [int(hex_color[i : i + 2], 16) for i in (0, 2, 4)]
    rgb = [int(c + (255 - c) * amount) for c in rgb]
    return "#" + "".join(f"{c:02x}" for c in rgb)


def darken(hex_color: str, amount: float = 0.2) -> str:
    hex_color = hex_color.lstrip("#")
    rgb = [int(hex_color[i : i + 2], 16) for i in (0, 2, 4)]
    rgb = [int(c * (1 - amount)) for c in rgb]
    return "#" + "".join(f"{c:02x}" for c in rgb)
