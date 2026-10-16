"""Legacy helper; prefer `pymagical dashboard --input_dir RESULTS_DIR`."""

import argparse
from pathlib import Path
import sys

# Support direct invocation from a source checkout or an installed package.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from pymagical.dashboard import add_arguments, launch_dashboard


def main():
    parser = argparse.ArgumentParser(description="Build the MAGICAL circuit dashboard")
    add_arguments(parser)
    launch_dashboard(parser.parse_args(), parser)


if __name__ == "__main__":
    main()
