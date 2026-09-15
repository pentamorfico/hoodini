"""Minimal Chrome DevTools Protocol (CDP) client for driving lightpanda.

Lightpanda (https://lightpanda.io/) is a lightweight, dependency-free headless
browser distributed via bioconda. It exposes a standard CDP WebSocket server,
so a small hand-rolled client built on ``websocket-client`` is enough to
navigate pages, fill in forms and read back page state -- no Playwright, no
Firefox, no GTK/X11 stack required.

This module only implements the tiny slice of CDP hoodini needs (navigate a
page, evaluate JS, wait for the load event) and is intentionally minimal.
"""

from __future__ import annotations

import contextlib
import json
import shutil
import subprocess
import time
import urllib.request
from collections.abc import Sequence

import websocket

LIGHTPANDA_WS_URL = "ws://127.0.0.1:9222/"
LIGHTPANDA_HTTP_URL = "http://127.0.0.1:9222/"
DEFAULT_COMMAND_TIMEOUT = 60


def find_lightpanda_binary() -> str | None:
    """Locate the ``lightpanda`` executable on PATH."""
    return shutil.which("lightpanda")


def lightpanda_server_alive(http_url: str = LIGHTPANDA_HTTP_URL) -> bool:
    """Check whether a lightpanda CDP server is already listening."""
    try:
        with urllib.request.urlopen(  # noqa: S310 -- fixed, hardcoded local CDP endpoint
            http_url + "json/version", timeout=3
        ) as resp:
            return resp.status == 200
    except Exception:
        return False


def start_lightpanda_server(lp_bin: str, wait_seconds: float = 15.0) -> bool:
    """Start ``lightpanda serve`` in the background if it's not already running.

    Returns True if a new server was started, False if one was already alive.
    """
    if lightpanda_server_alive():
        return False
    subprocess.Popen(
        [lp_bin, "serve"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True,
    )
    deadline = time.monotonic() + wait_seconds
    while time.monotonic() < deadline:
        if lightpanda_server_alive():
            return True
        time.sleep(0.5)
    raise RuntimeError(f"lightpanda serve did not come up on {LIGHTPANDA_HTTP_URL}")


class CDPError(RuntimeError):
    """Raised when a CDP command returns an error or a JS evaluation throws."""


class CDPSession:
    """A minimal synchronous CDP client over a single WebSocket connection."""

    def __init__(self, ws_url: str = LIGHTPANDA_WS_URL, timeout: float = DEFAULT_COMMAND_TIMEOUT):
        self.ws = websocket.create_connection(ws_url, timeout=timeout, suppress_origin=True)
        self._next_id = 1

    def send(
        self,
        method: str,
        params: dict | None = None,
        timeout: float = DEFAULT_COMMAND_TIMEOUT,
        session_id: str | None = None,
    ) -> dict:
        """Send a CDP command and block until its matching response arrives."""
        self.ws.settimeout(timeout)
        msg_id = self._next_id
        self._next_id += 1
        payload = {"id": msg_id, "method": method, "params": params or {}}
        if session_id:
            payload["sessionId"] = session_id
        self.ws.send(json.dumps(payload))
        while True:
            raw = self.ws.recv()
            msg = json.loads(raw)
            if msg.get("id") == msg_id:
                if "error" in msg:
                    raise CDPError(json.dumps(msg["error"]))
                return msg.get("result", {})
            # Ignore unrelated events/notifications while awaiting our reply.

    def wait_for_event(self, method: str, timeout: float = 30.0) -> bool:
        """Block until a CDP event of the given method name is received."""
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            try:
                self.ws.settimeout(1.0)
                raw = self.ws.recv()
                msg = json.loads(raw)
                if msg.get("method") == method:
                    return True
            except (TimeoutError, websocket.WebSocketTimeoutException):
                pass
        return False

    def open_page(self, timeout: float = 60.0) -> str:
        """Create an isolated browser context + target and return its session id."""
        browser_context = self.send("Target.createBrowserContext")
        target = self.send(
            "Target.createTarget",
            {"url": "about:blank", "browserContextId": browser_context["browserContextId"]},
        )
        attached = self.send(
            "Target.attachToTarget", {"targetId": target["targetId"], "flatten": True}
        )
        session_id = attached["sessionId"]
        self.send("Page.enable", {}, timeout, session_id)
        self.send("Runtime.enable", {}, timeout, session_id)
        return session_id

    def navigate(self, session_id: str, url: str, timeout: float = 120.0) -> None:
        """Navigate the page and wait for the load event to fire."""
        self.send("Page.navigate", {"url": url}, timeout, session_id)
        self.wait_for_event("Page.loadEventFired", timeout)

    def evaluate(self, session_id: str, expression: str, timeout: float = 30.0):
        """Evaluate JS in the page and return the resulting value."""
        result = self.send(
            "Runtime.evaluate",
            {"expression": expression, "awaitPromise": True, "returnByValue": True},
            timeout,
            session_id,
        )
        if result.get("exceptionDetails"):
            raise CDPError(json.dumps(result["exceptionDetails"]))
        return result.get("result", {}).get("value")

    def close(self) -> None:
        with contextlib.suppress(Exception):
            self.ws.close()

    def __enter__(self) -> CDPSession:
        return self

    def __exit__(self, *exc_info: Sequence) -> None:
        self.close()


__all__ = [
    "CDPError",
    "CDPSession",
    "LIGHTPANDA_HTTP_URL",
    "LIGHTPANDA_WS_URL",
    "find_lightpanda_binary",
    "lightpanda_server_alive",
    "start_lightpanda_server",
]
