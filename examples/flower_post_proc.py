#!/usr/bin/env python3
import argparse
from pathlib import Path
import sys
import inspect
# from convergence_study import plot_convergence_study_func
# from plot_flower import plot_all_films_func

import convergence_study
import plot_flower

# ----------------------------
# CLI arguments
# ----------------------------
parser = argparse.ArgumentParser(description="Run Flower post-processing on a case")
parser.add_argument(
    "yaml_file",
    type=str,
    help="Path to the YAML case file (e.g., ../Flower.jl/examples/one_fluid_hysing_coupled2.yml)"
)
parser.add_argument(
    "--func",
    type=str,
    default="plot_convergence_study_func",
    help="Function name to run (default: plot_convergence_study_func)"
)
# Add H5 files argument
parser.add_argument(
    "--h5",
    nargs="+",          # one or more .h5 files
    help="List of .h5 files to process"
)
parser.add_argument(
    "--skip-existing",
    action="store_true",
    help="Skip processing if the output figure already exists"
)
parser.add_argument(
    "--skiph5",
    nargs="+",          # one or more .h5 files
    help="Skip"
)
parser.add_argument(
    "--name",
    nargs="+",          # one or more figure names
    help="List of figures to process (for example plot only pressure_zoom defined in .yml in figures section and skip all others)"
)
parser.add_argument(
    "--dark",
    action="store_true",
    help=""
)
args = parser.parse_args()

# ----------------------------
# Normalize paths
# ----------------------------
script_dir = Path(__file__).resolve().parent     # directory of flower_post_proc.py
yaml_path = (Path(args.yaml_file)).expanduser().resolve()

if not yaml_path.exists():
    raise FileNotFoundError(f"YAML file not found: {yaml_path}")

# ----------------------------
# Import and call dynamically
# ----------------------------
# import convergence_study

if args.func in dir(convergence_study):
    func_to_run = getattr(convergence_study, args.func)
elif args.func in dir(plot_flower):
    func_to_run = getattr(plot_flower, args.func)
else:
    raise ValueError(f"No such function {args.func}")

# func_to_run = getattr(convergence_study, args.func)

sig = inspect.signature(func_to_run)

# Pass both yaml and args if function accepts them
if len(sig.parameters) == 2:
    func_to_run(yaml_path, args)
elif len(sig.parameters) == 1:
    func_to_run(yaml_path)
else:
    func_to_run()

print(f"✅ Ran {args.func} on {yaml_path}")
