#!/usr/bin/env python3
"""Generate bioconda meta.yaml from pyproject.toml.

Usage:
    python scripts/sync_from_pixi.py              # Generate meta.yaml
    python scripts/sync_from_pixi.py --check      # Check all deps
    python scripts/sync_from_pixi.py --check orfipy taxoniq  # Check specific deps
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import stat
import subprocess
import sys
import tarfile
import urllib.error
import urllib.request
from pathlib import Path

try:
    import tomllib
except ImportError:
    import tomli as tomllib  # type: ignore[import-not-found]

ROOT = Path(__file__).parent.parent
PYPROJECT_PATH = ROOT / "pyproject.toml"
META_PATH = ROOT / "recipes" / "meta.yaml"
MICROMAMBA_DIR = ROOT / ".micromamba"
MICROMAMBA_BIN = MICROMAMBA_DIR / "micromamba"

GITHUB_ORG = "pentamorfico"

# Build-only deps, not runtime
EXCLUDE_DEPS = {"pip", "setuptools", "wheel", "hoodini"}

# Platforms to check will be loaded from pyproject.toml

# Cache for conda package lookups
_conda_cache: dict[tuple[str, str, str], bool] = {}


def get_micromamba_url() -> str:
    """Get the correct micromamba download URL for this platform."""
    system = platform.system().lower()
    machine = platform.machine().lower()
    
    if system == "darwin":
        if machine in ("arm64", "aarch64"):
            return "https://micro.mamba.pm/api/micromamba/osx-arm64/latest"
        return "https://micro.mamba.pm/api/micromamba/osx-64/latest"
    elif system == "linux":
        if machine in ("arm64", "aarch64"):
            return "https://micro.mamba.pm/api/micromamba/linux-aarch64/latest"
        return "https://micro.mamba.pm/api/micromamba/linux-64/latest"
    else:
        raise RuntimeError(f"Unsupported platform: {system}-{machine}")


def ensure_micromamba() -> Path:
    """Download micromamba if not present, return path to binary."""
    if MICROMAMBA_BIN.exists():
        return MICROMAMBA_BIN
    
    print("Downloading micromamba...")
    MICROMAMBA_DIR.mkdir(parents=True, exist_ok=True)
    
    url = get_micromamba_url()
    tar_path = MICROMAMBA_DIR / "micromamba.tar.bz2"
    
    # Download
    with urllib.request.urlopen(url, timeout=60) as response:  # noqa: S310
        tar_path.write_bytes(response.read())
    
    # Extract
    with tarfile.open(tar_path, "r:bz2") as tar:
        for member in tar.getmembers():
            if member.name.endswith("micromamba"):
                member.name = "micromamba"
                tar.extract(member, MICROMAMBA_DIR, filter="data")
                break
    
    # Make executable
    MICROMAMBA_BIN.chmod(MICROMAMBA_BIN.stat().st_mode | stat.S_IEXEC)
    tar_path.unlink()
    
    print(f"✓ Installed micromamba to {MICROMAMBA_BIN}")
    return MICROMAMBA_BIN


def check_conda_available(
    package: str, 
    version_spec: str, 
    micromamba: Path,
    platforms: list[str]
) -> dict[str, bool]:
    """Check if a package (with version) is available on conda for each platform."""
    pkg_normalized = package.lower().replace("_", "-")
    
    # Build the search spec
    if version_spec:
        search_spec = f"{pkg_normalized}{version_spec}"
    else:
        search_spec = pkg_normalized
    
    results = {}
    for plat in platforms:
        cache_key = (pkg_normalized, version_spec, plat)
        if cache_key in _conda_cache:
            results[plat] = _conda_cache[cache_key]
            continue
        
        try:
            result = subprocess.run(
                [
                    str(micromamba), "search", search_spec,
                    "-c", "conda-forge", "-c", "bioconda",
                    "--platform", plat,
                    "--json", "-q"
                ],
                capture_output=True,
                text=True,
                timeout=30,
                env={**os.environ, "MAMBA_NO_BANNER": "1"}
            )
            
            # Parse JSON output
            available = False
            if result.returncode == 0 and result.stdout.strip():
                try:
                    data = json.loads(result.stdout)
                    pkgs = data.get("result", {}).get("pkgs", [])
                    available = len(pkgs) > 0
                except json.JSONDecodeError:
                    available = False
            
            _conda_cache[cache_key] = available
            results[plat] = available
        except (subprocess.TimeoutExpired, subprocess.SubprocessError):
            _conda_cache[cache_key] = False
            results[plat] = False
    
    return results


def parse_pep508(dep: str) -> tuple[str, str]:
    """Parse 'package>=1.0,<2' into ('package', '>=1.0,<2')."""
    match = re.match(r"^([a-zA-Z0-9_-]+)\s*(.*)", dep.strip())
    if match:
        return match.group(1), match.group(2)
    return dep.strip(), ""


def load_pyproject() -> dict:
    with open(PYPROJECT_PATH, "rb") as f:
        return tomllib.load(f)


def run_checks(
    all_deps: dict[str, str],
    linux_deps: dict[str, str],
    micromamba: Path,
    platforms: list[str],
    filter_pkgs: list[str] | None = None
) -> list[tuple[str, str, list[str]]]:
    """Run conda availability checks. Returns list of issues."""
    issues = []
    
    # Filter deps if specific packages requested
    if filter_pkgs:
        check_deps = {k: v for k, v in all_deps.items() if k.lower() in [p.lower() for p in filter_pkgs]}
        check_linux = {k: v for k, v in linux_deps.items() if k.lower() in [p.lower() for p in filter_pkgs]}
    else:
        check_deps = all_deps
        check_linux = linux_deps
    
    if not check_deps and not check_linux:
        print("No matching packages to check.")
        return issues
    
    # Check universal deps
    if check_deps:
        print(f"Checking conda availability across {platforms}...\n")
        header = f"{'Package':<25} {'Version':<20}"
        for plat in platforms:
            header += f" {plat:<14}"
        print(header)
        print("-" * (45 + 14 * len(platforms)))
        
        for pkg in sorted(check_deps.keys()):
            if pkg == "python":
                continue
            
            spec = check_deps[pkg]
            availability = check_conda_available(pkg, spec, micromamba, platforms)
            
            statuses = []
            missing = []
            for plat in platforms:
                if availability.get(plat, False):
                    statuses.append("✓")
                else:
                    statuses.append("✗")
                    missing.append(plat)
            
            spec_display = spec if spec else "*"
            line = f"{pkg:<25} {spec_display:<20}"
            for s in statuses:
                line += f" {s:<14}"
            print(line)
            
            if missing:
                issues.append((pkg, spec, missing))
    
    # Check Linux-specific deps
    if check_linux:
        print("\nLinux-specific dependencies:")
        linux_platforms = [p for p in platforms if p.startswith("linux")]
        if not linux_platforms:
            print("(No Linux platforms in checked set)")
        else:
            header = f"{'Package':<25} {'Version':<20}"
            for plat in linux_platforms:
                header += f" {plat:<14}"
            print(header)
            print("-" * (45 + 14 * len(linux_platforms)))
            
            for pkg in sorted(check_linux.keys()):
                spec = check_linux[pkg]
                availability = check_conda_available(pkg, spec, micromamba, platforms=linux_platforms)
                
                statuses = []
                missing = []
                for plat in linux_platforms:
                    if availability.get(plat, False):
                        statuses.append("✓")
                    else:
                        statuses.append("✗")
                        missing.append(plat)
                
                spec_display = spec if spec else "*"
                line = f"{pkg:<25} {spec_display:<20}"
                for s in statuses:
                    line += f" {s:<14}"
                print(line)
            
            if missing:
                issues.append((pkg + " [linux]", spec, missing))
    
    return issues


def generate_meta_yaml(config: dict) -> tuple[str, list[str], dict[str, str], dict[str, str]]:
    """Generate bioconda meta.yaml. Returns (yaml, git_deps, all_deps, linux_deps)."""
    project = config.get("project", {})
    pixi = config.get("tool", {}).get("pixi", {})

    name = project.get("name", "hoodini")
    version = project.get("version", "0.0.0")
    description = project.get("description", "")

    # Collect ALL dependencies
    all_deps: dict[str, str] = {}
    all_deps["python"] = ">=3.10"

    # From [project.dependencies]
    for dep in project.get("dependencies", []):
        pkg, spec = parse_pep508(dep)
        if pkg.lower() in EXCLUDE_DEPS:
            continue
        all_deps[pkg] = spec

    # From [tool.pixi.dependencies]
    for pkg, spec in pixi.get("dependencies", {}).items():
        if pkg.lower() in EXCLUDE_DEPS:
            continue
        if isinstance(spec, str):
            all_deps[pkg] = "" if spec == "*" else spec

    # Collect Linux-specific deps
    linux_deps: dict[str, str] = {}
    linux_target = pixi.get("target", {}).get("linux-64", {}).get("dependencies", {})
    for pkg, spec in linux_target.items():
        if isinstance(spec, str):
            linux_deps[pkg] = "" if spec == "*" else spec

    # Collect git dependencies
    git_deps = []
    for pkg, spec in pixi.get("pypi-dependencies", {}).items():
        if isinstance(spec, dict) and "git" in spec:
            git_deps.append(pkg)

    # Format deps for YAML
    run_deps = []
    for pkg in sorted(all_deps.keys()):
        spec = all_deps[pkg]
        formatted = f"{pkg} {spec}".strip() if spec else pkg
        run_deps.append(formatted)

    linux_deps_formatted = []
    for pkg in sorted(linux_deps.keys()):
        spec = linux_deps[pkg]
        formatted = f"{pkg} {spec}".strip() if spec else pkg
        linux_deps_formatted.append(f"{formatted}  # [linux]")

    run_deps_yaml = "\n    - ".join(run_deps)
    linux_deps_yaml = "\n    - ".join(linux_deps_formatted)

    # Get entry points
    scripts = project.get("scripts", {})
    entry_points_yaml = "\n    - ".join(
        f"{cmd} = {target}" for cmd, target in scripts.items()
    )

    meta = f'''# AUTO-GENERATED from pyproject.toml
{{% set name = "{name}" %}}
{{% set version = "{version}" %}}

package:
  name: {{{{ name|lower }}}}
  version: {{{{ version }}}}

source:
  url: https://github.com/{GITHUB_ORG}/{{{{ name }}}}/archive/refs/tags/v{{{{ version }}}}.tar.gz
  sha256: REPLACE_WITH_SHA256

build:
  number: 0
  noarch: python
  script: {{{{ PYTHON }}}} -m pip install . -vvv --no-deps --no-build-isolation
  entry_points:
    - {entry_points_yaml}

requirements:
  host:
    - python >=3.10
    - pip
    - setuptools >=68
  run:
    - {run_deps_yaml}'''

    if linux_deps_yaml:
        meta += f'''
    # Linux-specific
    - {linux_deps_yaml}'''

    meta += f'''

test:
  imports:
    - {name}
  commands:
    - {name} --help

about:
  home: https://github.com/{GITHUB_ORG}/{name}
  license: MIT
  license_family: MIT
  license_file: LICENSE
  summary: "{description}"
  dev_url: https://github.com/{GITHUB_ORG}/{name}

extra:
  recipe-maintainers:
    - pentamorfico
'''
    return meta, git_deps, all_deps, linux_deps


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate bioconda meta.yaml from pyproject.toml"
    )
    parser.add_argument(
        "--check",
        nargs="*",
        metavar="PKG",
        help="Check conda availability. No args = check all, or specify packages to check."
    )
    args = parser.parse_args()

    if not PYPROJECT_PATH.exists():
        print(f"Error: {PYPROJECT_PATH} not found", file=sys.stderr)
        return 1

    print(f"\nReading {PYPROJECT_PATH}...\n")
    config = load_pyproject()

    project = config.get("project", {})
    version = project.get("version", "0.0.0")

    # Generate meta.yaml
    meta_content, git_deps, all_deps, linux_deps = generate_meta_yaml(config)
    
    # Run checks if requested
    issues = []
    if args.check is not None:
        micromamba = ensure_micromamba()
        
        # Get platforms from config
        pixi = config.get("tool", {}).get("pixi", {})
        platforms = pixi.get("workspace", {}).get("platforms", [])
        if not platforms:
            platforms = ["linux-64", "linux-aarch64", "osx-64", "osx-arm64"]
            print("⚠ No platforms found in pyproject.toml, filtering to defaults.")

        filter_pkgs = args.check if args.check else None
        issues = run_checks(all_deps, linux_deps, micromamba, platforms, filter_pkgs)
        print()

    # Write meta.yaml
    META_PATH.parent.mkdir(parents=True, exist_ok=True)
    META_PATH.write_text(meta_content)
    print(f"✓ Generated {META_PATH}")

    # Show git deps warning only when not filtering specific packages
    show_git_warning = args.check is None or (args.check is not None and not args.check)
    if git_deps and show_git_warning:
        print(f"\n⚠ Git dependencies need conda OR bioconda recipes:")
        for dep in git_deps:
            print(f"   - {dep}")

    if issues:
        print(f"\n⚠ Packages missing on some platforms:")
        for pkg, spec, missing in issues:
            spec_str = f" {spec}" if spec else ""
            print(f"   - {pkg}{spec_str}: missing on {', '.join(missing)}")

    print(f"\nNext steps (if all dependencies are OK):")
    print(f"   1. Create a GitHub release/tag: git tag v{version} && git push --tags")
    print(f"   2. Get sha256: curl -sL <release_url>.tar.gz | shasum -a 256")
    print(f"   3. Replace REPLACE_WITH_SHA256 in meta.yaml")

    return 0


if __name__ == "__main__":
    sys.exit(main())
