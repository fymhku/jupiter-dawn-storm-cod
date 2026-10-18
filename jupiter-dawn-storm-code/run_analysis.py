"""Portable entry point for the 2017 manuscript analysis and figure sources."""
from pathlib import Path
import argparse
import os
import shutil
import subprocess
import sys
ROOT=Path(__file__).resolve().parent
TASKS={
 "check":"Verify input checksums and array contracts",
 "figures":"Replot Fig3, Fig4, S1 and S3 from included reference data",
 "psd":"Recompute calibrated PSD and numerical validation for both events",
 "cwt":"Recompute CWT amplitudes separately from the archived arrays",
 "bandpass":"Recompute and plot S1",
 "theory":"Run the original cubic and eighth-degree root scripts",
 "all":"Recompute PSD, CWT, bandpass, figures and original theory scripts",
}
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task",choices=list(TASKS),nargs="?",default="check")
    parser.add_argument("--list",action="store_true")
    parser.add_argument("--output",type=Path,default=ROOT/"output")
    parser.add_argument("--mag-dir",type=Path,help="Optional prepared daily JSS/RTP files for PSD; otherwise use included field windows")
    parser.add_argument("--save-coefficients",action="store_true")
    args=parser.parse_args()
    if args.list:
        for task,description in TASKS.items(): print(f"{task}: {description}")
        return
    output=args.output.resolve()
    if output==ROOT or output==ROOT/"data":
        parser.error("Choose a generated-output directory distinct from the package and its data directory")
    output.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy()
    env.update(JUPITER_RUN_ROOT=str(output),MPLBACKEND="Agg",PYTHONUTF8="1",PYTHONDONTWRITEBYTECODE="1")
    def run(script,options=(),cwd=None):
        print("Running "+script,flush=True)
        subprocess.run([sys.executable,"-X","utf8",str(ROOT/script),*map(str,options)],cwd=cwd or ROOT,env=env,check=True)
    if args.task=="check":
        run("code/check_package.py")
        return
    # Working copies can be regenerated, leaving the archived release inputs intact.
    target=output/"data"
    target.mkdir(parents=True,exist_ok=True)
    for path in (ROOT/"data").iterdir():
        dest=target/path.name
        if path.is_file() and not dest.exists(): shutil.copyfile(path,dest)
    if args.task in ["psd","all"]:
        options=[]
        if args.mag_dir: options += ["--mag-dir",args.mag_dir.resolve()]
        if args.save_coefficients: options += ["--save-coefficients"]
        run("code/psd_calibration.py",options)
    if args.task in ["cwt","all"]: run("code/recompute_cwt.py")
    if args.task in ["bandpass","figures","all"]: run("code/plot_bandpass.py")
    if args.task in ["figures","all"]:
        run("code/plot_fig3.py")
        run("code/plot_psd_figures.py")
    if args.task in ["theory","all"]:
        theory=output/"theory"
        theory.mkdir(exist_ok=True)
        for script in ["solve.py","solve_8.py"]: run("original/theory/"+script,cwd=theory)
    print(f"Completed. Outputs: {output}")
if __name__=="__main__": main()
