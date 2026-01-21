"""Memory tracking utilities for hoodini."""

from __future__ import annotations

import os
import threading
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Generator

try:
    import psutil

    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False

from hoodini.utils.logging_utils import console, debug, info


@dataclass
class MemoryStats:
    """Container for memory statistics."""

    current_mb: float = 0.0
    peak_mb: float = 0.0
    stage_stats: dict[str, float] = field(default_factory=dict)


def get_current_memory_mb() -> float:
    """Get current memory usage of the process in MB."""
    if not PSUTIL_AVAILABLE:
        return 0.0
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / (1024 * 1024)


def get_memory_info() -> dict[str, float]:
    """Get detailed memory information."""
    if not PSUTIL_AVAILABLE:
        return {"rss_mb": 0.0, "vms_mb": 0.0, "percent": 0.0}
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return {
        "rss_mb": mem_info.rss / (1024 * 1024),
        "vms_mb": mem_info.vms / (1024 * 1024),
        "percent": process.memory_percent(),
    }


def format_memory(mb: float) -> str:
    """Format memory value with appropriate units."""
    if mb >= 1024:
        return f"{mb / 1024:.2f} GB"
    return f"{mb:.1f} MB"


class MemoryTracker:
    """Track memory usage across pipeline stages."""

    def __init__(self, enabled: bool = True, poll_interval: float = 0.5):
        self.enabled = enabled and PSUTIL_AVAILABLE
        self.poll_interval = poll_interval
        self.stats = MemoryStats()
        self._monitor_thread: threading.Thread | None = None
        self._stop_event = threading.Event()
        self._lock = threading.Lock()

    def start_monitoring(self) -> None:
        """Start background memory monitoring thread."""
        if not self.enabled:
            return
        self._stop_event.clear()
        self._monitor_thread = threading.Thread(target=self._monitor_loop, daemon=True)
        self._monitor_thread.start()
        debug("Memory monitoring started")

    def stop_monitoring(self) -> None:
        """Stop background memory monitoring thread."""
        if not self.enabled or self._monitor_thread is None:
            return
        self._stop_event.set()
        self._monitor_thread.join(timeout=2.0)
        self._monitor_thread = None
        debug("Memory monitoring stopped")

    def _monitor_loop(self) -> None:
        """Background loop to track peak memory."""
        while not self._stop_event.is_set():
            current = get_current_memory_mb()
            with self._lock:
                self.stats.current_mb = current
                self.stats.peak_mb = max(self.stats.peak_mb, current)
            time.sleep(self.poll_interval)

    @contextmanager
    def track_stage(self, stage_name: str) -> Generator[None, None, None]:
        """Context manager to track memory for a pipeline stage."""
        if not self.enabled:
            yield
            return

        start_mem = get_current_memory_mb()
        peak_during_stage = start_mem

        # Local monitoring for this stage
        stage_stop = threading.Event()

        def stage_monitor():
            nonlocal peak_during_stage
            while not stage_stop.is_set():
                current = get_current_memory_mb()
                peak_during_stage = max(peak_during_stage, current)
                time.sleep(self.poll_interval)

        monitor = threading.Thread(target=stage_monitor, daemon=True)
        monitor.start()

        try:
            yield
        finally:
            stage_stop.set()
            monitor.join(timeout=1.0)
            end_mem = get_current_memory_mb()

            with self._lock:
                self.stats.stage_stats[stage_name] = peak_during_stage
                self.stats.peak_mb = max(self.stats.peak_mb, peak_during_stage)

            delta = end_mem - start_mem
            delta_str = f"+{delta:.1f}" if delta >= 0 else f"{delta:.1f}"
            debug(
                f"Stage '{stage_name}': peak={format_memory(peak_during_stage)}, "
                f"delta={delta_str} MB"
            )

    def get_peak_memory(self) -> float:
        """Get peak memory usage in MB."""
        with self._lock:
            return self.stats.peak_mb

    def get_stage_stats(self) -> dict[str, float]:
        """Get memory stats per stage."""
        with self._lock:
            return dict(self.stats.stage_stats)

    def print_summary(self) -> None:
        """Print a summary of memory usage."""
        if not self.enabled:
            info("Memory tracking disabled (psutil not available)")
            return

        peak = self.get_peak_memory()
        stage_stats = self.get_stage_stats()

        console.print()
        console.print("[bold cyan]📊 Memory Usage Summary[/bold cyan]")
        console.print(f"  Peak memory: [bold]{format_memory(peak)}[/bold]")

        if stage_stats:
            console.print("  Per-stage peak memory:")
            # Sort by memory usage (descending)
            sorted_stages = sorted(stage_stats.items(), key=lambda x: x[1], reverse=True)
            for stage, mem in sorted_stages:
                console.print(f"    • {stage}: {format_memory(mem)}")
        console.print()


# Global tracker instance
_global_tracker: MemoryTracker | None = None


def get_tracker() -> MemoryTracker:
    """Get or create the global memory tracker."""
    global _global_tracker
    if _global_tracker is None:
        _global_tracker = MemoryTracker()
    return _global_tracker


def reset_tracker(enabled: bool = True) -> MemoryTracker:
    """Reset and return a new global memory tracker."""
    global _global_tracker
    _global_tracker = MemoryTracker(enabled=enabled)
    return _global_tracker
