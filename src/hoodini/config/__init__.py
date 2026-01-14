from importlib.resources import files

import tomli

from hoodini.config.settings import RuntimeConfig, build_runtime_config


def load_default_config() -> dict:
    defaults_path = files("hoodini.config").joinpath("defaults.toml")
    with open(defaults_path, "rb") as f:
        grouped = tomli.load(f)

    flat = {}
    for section, values in grouped.items():
        if isinstance(values, dict):
            flat.update(values)
        else:
            flat[section] = values
    return flat


__all__ = ["RuntimeConfig", "build_runtime_config", "load_default_config"]
