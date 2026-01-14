"""Deprecated Config schema.

Use `hoodini.config.settings.RuntimeConfig` for all new code. This alias exists
to avoid breaking imports while the codebase migrates to the new config layer.
"""

from __future__ import annotations

from hoodini.config.settings import RuntimeConfig as Config

__all__ = ["Config"]
