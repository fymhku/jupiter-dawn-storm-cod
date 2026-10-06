"""Run the supported workflows from any working directory."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
TASKS = {
    "field": "juno_ulf_analysis/run_field_plot.py",
    "cwt": "juno_ulf_analysis/run_cwt.py",
    "psd": "juno_ulf_analysis/run_psd.py",
    "bandpass": "juno_ulf_analysis/run_bandpass.py",
    "quiet": "supplementary/make_quiet_figs.py",
    "quiet-panels": "supplementary/stack_0510.py",
    "theory": "theory/solve.py",
    "theory-8": "theory/solve_8.py",
}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", choices=[*TASKS, "all"], nargs="?", default="all")
    parser.add_argument("--list", action="store_true", help="list tasks without importing scientific libraries")
    args = parser.parse_args()
    if args.list:
        for task, file in TASKS.items():
            print(f"{task}: {file}")
        return
    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")
    env.setdefault("PYTHONUTF8", "1")
    for task in (TASKS if args.task == "all" else [args.task]):
        cwd = ROOT
        if task.startswith("theory"):
            cwd = Path(env.get("JUNO_ULF_OUT", str(ROOT / "output"))).resolve() / "theory"
            cwd.mkdir(parents=True, exist_ok=True)
        print(f"Running {task}", flush=True)
        subprocess.run([sys.executable, str(ROOT / TASKS[task])],
                       cwd=cwd, env=env, check=True)

if __name__ == "__main__":
    main()
