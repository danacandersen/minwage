#!/usr/bin/env python3
"""
run_analysis.py
Planter Minimum Wage Project

Master runner: calls Stata on each do-file in sequence (or a specific step).

Usage:
    python run_analysis.py              # run all steps 00-03
    python run_analysis.py --step 00    # run only 00_data_prep.do
    python run_analysis.py --step 02    # run only 02_main_analysis.do
    python run_analysis.py --from 01    # run steps 01 through 03

Stata executable is auto-detected from common macOS paths.
Override with --stata /path/to/stata.

Log files are saved by the do-files themselves to logs/*.log.
"""

import argparse
import subprocess
import sys
import os
import time
from datetime import datetime
from pathlib import Path


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Do-files to run, in order. Keys are step numbers (as strings).
DO_FILES = {
    "00": "00_data_prep.do",
    "01": "01_data_checks.do",
    "02": "02_main_analysis.do",
    "03": "03_heterogeneity.do",
}

# Common Stata executable paths on macOS.
# Add Windows/Linux paths here if needed.
STATA_CANDIDATES = [
    # Stata 19 MP (primary — user's version)
    "/Applications/Stata19/StataMP.app/Contents/MacOS/StataMP",
    "/Applications/Stata/StataMP 19.app/Contents/MacOS/StataMP",
    "/Applications/StataMP 19.app/Contents/MacOS/StataMP",
    # Stata 19 generic location (some installers use this)
    "/Applications/Stata/StataMP.app/Contents/MacOS/StataMP",
    "/Applications/StataMP.app/Contents/MacOS/StataMP",
    # Stata 19 SE fallback
    "/Applications/Stata19/StataSE.app/Contents/MacOS/StataSE",
    "/Applications/Stata/StataSE.app/Contents/MacOS/StataSE",
    "/Applications/StataSE.app/Contents/MacOS/StataSE",
    # Stata 18
    "/Applications/Stata18/StataMP.app/Contents/MacOS/StataMP",
    "/Applications/Stata18/StataSE.app/Contents/MacOS/StataSE",
    # Stata 17
    "/Applications/Stata17/StataMP.app/Contents/MacOS/StataMP",
    "/Applications/Stata17/StataSE.app/Contents/MacOS/StataSE",
    # Generic (no version in path)
    "/Applications/Stata/Stata.app/Contents/MacOS/Stata",
    "/Applications/Stata.app/Contents/MacOS/Stata",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_stata(override=None) -> str:
    """Return path to Stata executable, or raise if not found."""
    if override:
        if os.path.isfile(override):
            return override
        raise FileNotFoundError(f"Stata not found at specified path: {override}")

    for path in STATA_CANDIDATES:
        if os.path.isfile(path):
            return path

    raise FileNotFoundError(
        "Stata executable not found. Checked:\n"
        + "\n".join(f"  {p}" for p in STATA_CANDIDATES)
        + "\n\nUse --stata /path/to/stata to specify manually."
    )


def ts() -> str:
    """Return timestamp string for console output."""
    return datetime.now().strftime("%H:%M:%S")


def run_do_file(stata_exe: str, do_file: str, project_dir: Path) -> int:
    """
    Run a single Stata do-file in batch mode.
    Returns the Stata process return code (0 = success).
    """
    do_path = project_dir / do_file

    if not do_path.exists():
        print(f"  [ERROR] Do-file not found: {do_path}")
        return 1

    # Pass only the filename (not full path) so spaces in the project
    # directory don't confuse Stata's batch-mode argument parser.
    # This works because cwd is already set to project_dir below.
    cmd = [stata_exe, "-b", "do", do_file]

    print(f"  Command: {' '.join(cmd)}")
    print(f"  Working dir: {project_dir}")

    try:
        result = subprocess.run(
            cmd,
            cwd=str(project_dir),
            capture_output=True,
            text=True,
            timeout=3600,  # 1-hour timeout per do-file
        )
    except subprocess.TimeoutExpired:
        print(f"  [ERROR] Timed out after 3600s.")
        return 1
    except Exception as e:
        print(f"  [ERROR] Failed to launch Stata: {e}")
        return 1

    # Stata batch mode writes output to a .log file with same stem as do-file;
    # also capture stdout/stderr for diagnostics.
    if result.stdout:
        print("  [stdout]", result.stdout[:500])
    if result.stderr:
        print("  [stderr]", result.stderr[:500])

    return result.returncode


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Run Planter MW Analysis do-files via Stata."
    )
    parser.add_argument(
        "--step",
        metavar="NN",
        help="Run only this step (e.g. 00, 01, 02, 03).",
    )
    parser.add_argument(
        "--from",
        dest="from_step",
        metavar="NN",
        help="Run all steps starting from this step.",
    )
    parser.add_argument(
        "--stata",
        metavar="PATH",
        help="Override Stata executable path.",
    )
    parser.add_argument(
        "--dir",
        metavar="PATH",
        default=None,
        help="Project directory (default: directory containing this script).",
    )
    args = parser.parse_args()

    # Project directory
    if args.dir:
        project_dir = Path(args.dir).resolve()
    else:
        project_dir = Path(__file__).parent.resolve()

    print(f"[{ts()}] Project dir: {project_dir}")

    # Detect Stata
    try:
        stata_exe = find_stata(args.stata)
        print(f"[{ts()}] Stata: {stata_exe}")
    except FileNotFoundError as e:
        print(f"[ERROR] {e}")
        sys.exit(1)

    # Determine which steps to run
    all_steps = sorted(DO_FILES.keys())

    if args.step:
        if args.step not in DO_FILES:
            print(f"[ERROR] Unknown step '{args.step}'. Valid: {', '.join(all_steps)}")
            sys.exit(1)
        steps = [args.step]
    elif args.from_step:
        if args.from_step not in DO_FILES:
            print(f"[ERROR] Unknown step '{args.from_step}'. Valid: {', '.join(all_steps)}")
            sys.exit(1)
        steps = [s for s in all_steps if s >= args.from_step]
    else:
        steps = all_steps

    print(f"[{ts()}] Steps to run: {steps}")
    print()

    # Create logs/ and output/ directories if missing
    (project_dir / "logs").mkdir(exist_ok=True)
    (project_dir / "output").mkdir(exist_ok=True)

    # Run each step
    results = {}
    overall_start = time.time()

    for step in steps:
        do_file = DO_FILES[step]
        print(f"[{ts()}] ── Step {step}: {do_file} ──────────────────────────────")

        step_start = time.time()
        rc = run_do_file(stata_exe, do_file, project_dir)
        elapsed = time.time() - step_start

        results[step] = rc
        status = "OK" if rc == 0 else f"FAILED (rc={rc})"
        print(f"[{ts()}] Step {step}: {status}  ({elapsed:.0f}s)")
        print()

        # Stop on failure (downstream steps depend on earlier ones)
        if rc != 0:
            print(f"[{ts()}] Stopping due to failure in step {step}.")
            break

    # Summary
    total = time.time() - overall_start
    print("=" * 60)
    print(f"[{ts()}] Summary ({total:.0f}s total)")
    for step, rc in results.items():
        status = "OK" if rc == 0 else f"FAILED (rc={rc})"
        print(f"  Step {step} ({DO_FILES[step]}): {status}")

    skipped = [s for s in steps if s not in results]
    for step in skipped:
        print(f"  Step {step} ({DO_FILES[step]}): SKIPPED (earlier failure)")

    print("=" * 60)
    print(f"[{ts()}] Logs: {project_dir / 'logs'}")
    print(f"[{ts()}] Figures: {project_dir / 'output'}")

    # Exit with failure code if any step failed
    if any(rc != 0 for rc in results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
