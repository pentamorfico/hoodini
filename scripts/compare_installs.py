import argparse
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Generator
from contextlib import contextmanager
from pathlib import Path

RED = "\033[91m"
GREEN = "\033[92m"
YELLOW = "\033[93m"
RESET = "\033[0m"

ROOT = Path(__file__).parent.parent
PIXI_LOCK = ROOT / "pixi.lock"

PLATFORM_MAP = {
    "linux-64": "x86_64-unknown-linux-gnu",
    "linux-aarch64": "aarch64-unknown-linux-gnu",
    "osx-64": "x86_64-apple-darwin",
    "osx-arm64": "aarch64-apple-darwin",
}


def get_current_platform() -> str:
    """Detect current platform string for pixi."""
    sys_name = platform.system().lower()
    machine = platform.machine().lower()

    if sys_name == "linux":
        return "linux-aarch64" if machine == "aarch64" else "linux-64"
    elif sys_name == "darwin":
        return "osx-arm64" if machine == "arm64" else "osx-64"
    return "unknown"


def parse_pixi_lock(lock_path: Path, target_platform: str) -> dict[str, str]:
    """Parse packages from a pixi.lock file for a specific platform."""
    if not lock_path.exists():
        return {}

    packages = {}
    current_section = None

    with open(lock_path) as f:
        for line in f:
            line = line.strip()
            if line.endswith(":") and not line.startswith("-"):
                current_section = line[:-1]
                continue

            if current_section == target_platform and line.startswith("- "):
                parts = line.split()
                if len(parts) < 3:
                    continue

                source_type = parts[1]
                url = parts[2]
                filename = url.split("/")[-1]

                for ext in [".conda", ".tar.bz2", ".whl", ".tar.gz"]:
                    if filename.endswith(ext):
                        filename = filename[: -len(ext)]
                        break

                if source_type == "pypi:":
                    match = re.match(
                        r"^([a-zA-Z0-9_\-]+?)-(\d+(?:\.\d+)*(?:[a-zA-Z0-9_\.]*)?)", filename
                    )
                    if match:
                        name = match.group(1).lower().replace("_", "-")
                        version = match.group(2)
                        packages[name] = version
                else:
                    match = re.search(
                        r"^(.+)-(\d+(?:\.\d+)*(?:[a-zA-Z0-9_\.]*)?)-([^-]+)$", filename
                    )
                    if match:
                        name = match.group(1).lower().replace("_", "-")
                        version = match.group(2)
                        packages[name] = version
    return packages


@contextmanager
def get_lock_file(fresh: bool) -> Generator[Path, None, None]:
    """Context manager yielding path to lock file (existing or temp fresh)."""
    if not fresh:
        if not PIXI_LOCK.exists():
            print(f"{RED}Error: pixi.lock not found{RESET}")
            sys.exit(1)
        print("Using existing pixi.lock")
        yield PIXI_LOCK
    else:
        print("Resolving fresh pixi environment (this may take a minute)...")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)

            # Copy necessary files
            for filename in ["pyproject.toml", "README.md", "README", "LICENSE", "LICENSE.md"]:
                src = ROOT / filename
                if src.exists():
                    shutil.copy(src, tmp_path / filename)

            # Copy src directory for editable install
            src_dir = ROOT / "src"
            if src_dir.exists():
                shutil.copytree(src_dir, tmp_path / "src")

            # Run pixi install
            try:
                subprocess.run(["pixi", "install"], cwd=tmp_path, check=True, capture_output=True)
            except subprocess.CalledProcessError as e:
                print(f"{RED}pixi install failed:{RESET}")
                print(e.stderr.decode() if e.stderr else "Unknown error")
                sys.exit(1)

            yield tmp_path / "pixi.lock"


def get_uv_resolved_packages(target_platform: str) -> dict[str, str]:
    """Resolve packages using uv pip compile for the target platform."""
    uv_platform = PLATFORM_MAP.get(target_platform)
    if not uv_platform:
        print(f"{RED}Unknown platform mapping for {target_platform}{RESET}")
        return {}

    cmd = ["uv", "pip", "compile", "pyproject.toml", "--python-platform", uv_platform, "--quiet"]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"{RED}uv resolution failed for {target_platform}:{RESET}")
        print(e.stderr)
        return {}

    packages = {}
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        match = re.match(r"^([a-zA-Z0-9_\-]+)==([^ ;]+)", line)
        if match:
            name = match.group(1).lower().replace("_", "-")
            version = match.group(2)
            packages[name] = version

    return packages


def compare_platform(platform: str, lock_path: Path) -> list[str]:
    """Compare packages for a single platform and return output lines."""
    lines = []
    lines.append(f"Checking {platform}...")

    pixi_pkgs = parse_pixi_lock(lock_path, platform)
    if not pixi_pkgs:
        lines.append(f"No packages found in pixi.lock for {platform}")
        return lines

    uv_pkgs = get_uv_resolved_packages(platform)

    all_pkgs = sorted(set(pixi_pkgs.keys()) | set(uv_pkgs.keys()))

    lines.append(f"\n{'Package':<30} {'Pixi Version':<20} {'UV Version':<20} {'Status':<20}")
    lines.append("-" * 90)

    mismatches = 0
    pixi_only = 0
    uv_only = 0
    mismatch_list = []

    for pkg in all_pkgs:
        p_ver = pixi_pkgs.get(pkg, "-")
        u_ver = uv_pkgs.get(pkg, "-")

        if p_ver == u_ver:
            status = "Match"
        elif p_ver == "-":
            status = "UV Only"
            uv_only += 1
        elif u_ver == "-":
            status = "Pixi Only"
            pixi_only += 1
        else:
            status = "Mismatch"
            mismatches += 1
            mismatch_list.append((pkg, p_ver, u_ver))

        # For logging, we strip colors. For stdout, we keep them if interactive?
        # Let's keep formatted string for return, and maybe colorize only when printing to term
        line = f"{pkg:<30} {p_ver:<20} {u_ver:<20} {status}"
        lines.append(line)

    lines.append("-" * 90)
    lines.append(f"Summary for {platform}:")
    lines.append(f"  Matches:    {len(all_pkgs) - mismatches - pixi_only - uv_only}")
    lines.append(f"  Mismatches: {mismatches}")
    lines.append(f"  Pixi Only:  {pixi_only}")
    lines.append(f"  UV Only:    {uv_only}")

    if mismatch_list:
        lines.append("\nMismatches:")
        lines.append(f"{'Package':<30} {'Pixi':<20} {'UV':<20}")
        lines.append("-" * 72)
        for pkg, p, u in mismatch_list:
            lines.append(f"{pkg:<30} {p:<20} {u:<20}")

    return lines


def main():
    parser = argparse.ArgumentParser(description="Compare Pixi lock vs UV resolution.")
    parser.add_argument(
        "--platform",
        default=get_current_platform(),
        choices=list(PLATFORM_MAP.keys()),
        help="Target platform to check",
    )
    parser.add_argument(
        "--all-platforms", action="store_true", help="Check all platforms defined in script"
    )
    parser.add_argument(
        "--use-lock-file",
        action="store_true",
        help="Use existing pixi.lock instead of fresh resolve",
    )
    args = parser.parse_args()

    platforms = list(PLATFORM_MAP.keys()) if args.all_platforms else [args.platform]

    with get_lock_file(fresh=not args.use_lock_file) as lock_path:
        for plat in platforms:
            output_lines = compare_platform(plat, lock_path)

            # Print to stdout (plain)
            print("\n".join(output_lines))
            print()

            # Write to log
            log_file = ROOT / f"install_check_{plat}.log"
            with open(log_file, "w") as f:
                f.write("\n".join(output_lines))
                f.write("\n")
            print(f"Report written to {log_file}")


if __name__ == "__main__":
    main()
