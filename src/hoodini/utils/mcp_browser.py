"""Minimal MCP-stdio client for driving headless browsers without ports.

Supports two drivers over the same JSON-RPC/MCP line protocol:

- ``obscura`` (preferred: ships on conda-forge), tools ``browser_navigate`` /
  ``browser_evaluate`` / ``browser_search``.
- ``lightpanda`` (bioconda package or the official pip wheel), tools ``goto`` /
  ``evaluate``.

Both are child processes talked to over stdin/stdout pipes: no CDP server and
no TCP port is involved, which is friendlier to HPC nodes with strict network
policies and avoids port collisions between concurrent hoodini runs.

This module implements the tiny slice hoodini needs (navigate, evaluate JS,
wait for a condition, search page text) on top of plain subprocess pipes — no
MCP SDK required.
"""

from __future__ import annotations

import contextlib
import json
import queue
import shutil
import subprocess
import threading
import time
from collections.abc import Sequence

DEFAULT_TIMEOUT = 120.0
MCP_PROTOCOL_VERSION = "2024-11-05"


def find_lightpanda_binary() -> str | None:
    """Locate a lightpanda binary: PATH first, then the pip wheel's bundled one."""
    found = shutil.which("lightpanda")
    if found:
        return found
    try:
        from lightpanda.client import find_binary as pip_find_binary

        return str(pip_find_binary())
    except ImportError:
        return None


def find_browser_binary(preferred: str = "obscura") -> tuple[str, str] | None:
    """Locate a supported headless browser binary.

    Returns ``(driver, path)`` where driver is ``"obscura"`` or
    ``"lightpanda"``, trying the preferred browser first and falling back to
    the other one. Obscura ships on conda-forge; lightpanda on bioconda and
    as the official ``lightpanda`` pip wheel.
    """
    order = [preferred, "lightpanda"] if preferred == "obscura" else ["lightpanda", "obscura"]
    for driver in order:
        found = shutil.which(driver)
        if found:
            return driver, found
    for driver in order:
        if driver == "lightpanda":
            try:
                from lightpanda.client import find_binary as pip_find_binary

                return "lightpanda", str(pip_find_binary())
            except ImportError:
                continue
    return None


class MCPBrowserError(RuntimeError):
    """Raised when a browser MCP tool call fails or times out."""


class MCPBrowserSession:
    """Drive a headless browser over MCP stdio (JSON-RPC line protocol, no TCP port)."""

    def __init__(
        self,
        binary: str | None = None,
        driver: str | None = None,
        timeout: float = DEFAULT_TIMEOUT,
    ):
        self.timeout = timeout
        self._next_id = 0
        resolved = find_browser_binary() if binary is None else ("obscura", binary)
        if binary is not None:
            driver = driver or "obscura"
        if resolved is None:
            raise MCPBrowserError("no obscura or lightpanda binary found on PATH")
        self.driver, self.binary = resolved
        self._proc = subprocess.Popen(  # noqa: S603 -- resolved binary path
            [self.binary, "mcp"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            bufsize=1,
        )
        self._queue: queue.Queue = queue.Queue()
        self._reader = threading.Thread(target=self._read_loop, daemon=True)
        self._reader.start()
        self._call(
            "initialize",
            {
                "protocolVersion": MCP_PROTOCOL_VERSION,
                "capabilities": {},
                "clientInfo": {"name": "hoodini", "version": "1.0"},
            },
        )
        self._notify("notifications/initialized")

    # low-level ----------------------------------------------------------

    def _write(self, msg: dict) -> None:
        assert self._proc.stdin is not None
        self._proc.stdin.write(json.dumps(msg) + "\n")
        self._proc.stdin.flush()

    def _read_loop(self) -> None:
        assert self._proc.stdout is not None
        for line in self._proc.stdout:
            self._queue.put(line)
        self._queue.put(None)  # EOF sentinel

    def _read(self, want_id: int, deadline: float) -> dict:
        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise MCPBrowserError(f"timed out waiting for MCP response id={want_id}")
            try:
                line = self._queue.get(timeout=min(remaining, 1.0))
            except queue.Empty:
                continue
            if line is None:
                raise MCPBrowserError("browser closed stdio unexpectedly")
            line = line.strip()
            if not line:
                continue
            msg = json.loads(line)
            if msg.get("id") == want_id:
                return msg

    def _call(self, method: str, params: dict | None = None, timeout: float | None = None) -> dict:
        self._next_id += 1
        msg_id = self._next_id
        msg: dict = {"jsonrpc": "2.0", "id": msg_id, "method": method}
        if params is not None:
            msg["params"] = params
        self._write(msg)
        return self._read(msg_id, time.monotonic() + (timeout or self.timeout))

    def _notify(self, method: str, params: dict | None = None) -> None:
        msg: dict = {"jsonrpc": "2.0", "method": method}
        if params is not None:
            msg["params"] = params
        self._write(msg)

    # tools ----------------------------------------------------------------

    @staticmethod
    def _tool_text(result: dict) -> str:
        parts = [
            block.get("text", "")
            for block in result.get("content", [])
            if block.get("type") == "text"
        ]
        return "\n".join(parts)

    def tool(self, name: str, arguments: dict | None = None, timeout: float | None = None) -> str:
        """Invoke a browser MCP tool and return its text content."""
        resp = self._call("tools/call", {"name": name, "arguments": arguments or {}}, timeout)
        result = resp.get("result", {})
        if result.get("isError"):
            raise MCPBrowserError(f"tool '{name}' failed: {self._tool_text(result)}")
        return self._tool_text(result)

    # driver-neutral high level ---------------------------------------------

    def goto(self, url: str, timeout: float | None = None) -> None:
        """Navigate the page and wait for it to load."""
        self.tool("browser_navigate" if self.driver == "obscura" else "goto", {"url": url}, timeout)

    def evaluate(self, script: str, timeout: float | None = None):
        """Evaluate JS in the page; return the parsed value when possible."""
        tool = "browser_evaluate" if self.driver == "obscura" else "evaluate"
        return self._parse_value(
            self.tool(
                tool, {"script" if self.driver == "lightpanda" else "expression": script}, timeout
            )
        )

    @staticmethod
    def _parse_value(text: str):
        """Parse a tool text payload into its JS value.

        lightpanda appends annotation lines like ``(Navigated to ...)`` after
        the value; try the whole text first, then shorter line prefixes, and
        fall back to the raw text.
        """
        text = text.strip()
        try:
            return json.loads(text)
        except ValueError:
            pass
        lines = text.splitlines()
        for end in range(1, len(lines) + 1):
            try:
                return json.loads("\n".join(lines[:end]))
            except ValueError:
                continue
        return text

    def search_text(self, query: str, limit: int = 3, timeout: float | None = None) -> str:
        """Search the rendered page text (obscura only)."""
        return self.tool("browser_search", {"query": query, "limit": limit}, timeout)

    @property
    def supports_search(self) -> bool:
        return self.driver == "obscura"

    def wait_for_script(self, script: str, timeout: float = 30.0, poll: float = 1.0) -> bool:
        """Poll a JS expression until it returns truthy (or timeout)."""
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            try:
                if self.evaluate(script):
                    return True
            except MCPBrowserError:
                pass
            time.sleep(poll)
        return False

    def close(self) -> None:
        with contextlib.suppress(Exception):
            self._proc.stdin.close()
        try:
            self._proc.wait(timeout=10)
        except Exception:
            with contextlib.suppress(Exception):
                self._proc.kill()

    def __enter__(self) -> MCPBrowserSession:
        return self

    def __exit__(self, *exc_info: Sequence) -> None:
        self.close()


__all__ = [
    "MCPBrowserError",
    "MCPBrowserSession",
    "find_browser_binary",
    "find_lightpanda_binary",
]
