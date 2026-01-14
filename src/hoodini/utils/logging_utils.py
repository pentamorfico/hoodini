"""Centralized logging utilities with Rich styling and log levels."""

from __future__ import annotations

import logging
from collections.abc import Callable
from datetime import datetime
from typing import Any

from rich.console import Console
from rich.panel import Panel

_QUIET = False
_DEBUG = False


def _make_console(quiet: bool, debug: bool) -> Console:
    # We add timestamps manually to messages to keep them visible during progress rendering.
    c = Console(
        quiet=quiet,
        log_time=False,
        log_path=False,
    )
    return c


console = _make_console(False, False)
logger = logging.getLogger("hoodini")


def configure_logging(*, quiet: bool = False, debug: bool = False) -> None:
    """Configure console/logging flags."""
    global _QUIET, _DEBUG, console
    _QUIET = quiet
    _DEBUG = debug
    console = _make_console(quiet, debug)
    level = logging.DEBUG if debug else logging.INFO
    logger.setLevel(level)


def info(message: str) -> None:
    if _QUIET:
        return
    console.print(f"{_ts()} [light_slate_grey]{message}[/light_slate_grey]")


def success(message: str) -> None:
    if _QUIET:
        return
    console.print(f"{_ts()} [green]✔ {message}[/green]")


def warn(message: str) -> None:
    if _QUIET:
        return
    console.print(f"{_ts()} [orange3]Warning:[/orange3] {message}")


def error(message: str) -> None:
    if _QUIET:
        return
    console.print(f"{_ts()} [bright_red]Error:[/bright_red] {message}")


def debug(message: str) -> None:
    if _QUIET or not _DEBUG:
        return
    console.print(f"{_ts()} [dim]{message}[/dim]")


def is_debug_enabled() -> bool:
    return _DEBUG and not _QUIET


def header(title: str, subtitle: str | None = None, border_style: str = "light_slate_grey") -> None:
    """Render a boxed header with optional subtitle."""
    if _QUIET:
        return
    text = f"[bold light_slate_grey]{title}[/bold light_slate_grey]"
    if subtitle:
        text += f"\n[dim]{subtitle}[/dim]"
    console.print(Panel.fit(text, border_style=border_style))


def stage_header(title: str, emoji: str = "") -> None:
    """Log a stage banner with timestamp."""
    header(f"{emoji} {title}" if emoji else title)


def stage_done(message: str) -> None:
    """Log a completion message with timestamp."""
    if _QUIET:
        return
    console.print(Panel.fit(f"[bold green]✔️  {message}[/bold green]", border_style="green"))
    console.print("[light_slate_grey]" + "─" * 80 + "[/light_slate_grey]")


def prompt(message: str, default: str | None = None) -> str:
    """Prompt user input aligned with log indentation."""
    if _QUIET:
        return default or ""
    suffix = f" ({default})" if default else ""
    prompt_text = f"{_ts()} {message}{suffix}: "
    return console.input(prompt_text).strip()


def run_with_spinner(
    title: str, func: Callable[..., Any], *args, spinner_name: str = "dots", **kwargs
):
    """
    Display a Rich spinner with `title` while running `func(*args, **kwargs)`.
    Returns whatever `func` returns.
    """
    if _QUIET:
        return func(*args, **kwargs)
    with console.status(f"[bold cyan]{title}[/bold cyan]", spinner=spinner_name):
        return func(*args, **kwargs)


def _ts() -> str:
    t = datetime.now().strftime("%H:%M:%S")
    return f"[grey53][[/grey53][light_slate_grey]{t}[/light_slate_grey][grey53]][/grey53]"
