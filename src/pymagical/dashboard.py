"""Build the bundled dashboard from a directory of MAGICAL results."""

import json
from importlib.resources import files
from pathlib import Path
import webbrowser


RESULT_COLUMNS = (
    "Gene_symbol", "Gene_chr", "Gene_TSS", "Peak_chr", "Peak_start",
    "Peak_end", "Peak_Gene_Prob",
)
DATA_SCRIPT = '<script src="./data.js"></script>'


def generate_dashboard(input_dir, output=None):
    """Save standalone HTML, returning its absolute path.

    Matrix and timing sidecars are ignored; installed package files stay intact.
    """
    directory = Path(input_dir).expanduser().resolve()
    if not directory.is_dir():
        raise ValueError(f"Input directory does not exist or is not a directory: {directory}")
    results = {}
    for path in sorted(directory.glob("*.txt")):
        if not path.is_file():
            continue
        with path.open(encoding="utf-8-sig") as handle:
            header = handle.readline()
            columns = header.rstrip("\r\n").split("\t")
            if tuple(columns[:7]) != RESULT_COLUMNS or len(columns) < 8:
                continue
            results[path.stem] = header + handle.read()
    if not results:
        raise ValueError(f"No MAGICAL result tables found in: {directory}")

    template = files("pymagical").joinpath("ui", "magical_explorer.html").read_text(encoding="utf-8")
    # Keep names/data containing '</script>' inside literal JS strings.
    payload = json.dumps(results, ensure_ascii=True).replace("<", "\\u003c")
    html = template.replace(DATA_SCRIPT, f"<script>window.MAGICAL_RESULTS = {payload};</script>")
    destination = Path(output).expanduser().resolve() if output else directory / "dashboard.html"
    destination.write_text(html, encoding="utf-8")
    return destination


def add_arguments(parser):
    """Share arguments between the package CLI and the legacy helper."""
    parser.add_argument(
        "--input_dir", "--input-dir", required=True,
        help="Directory containing pymagical output/result .txt files",
    )
    parser.add_argument("--output", help="Output HTML path (default: INPUT_DIR/dashboard.html)")
    parser.add_argument("--no-open", action="store_true", help="Save without opening a browser")


def launch_dashboard(args, parser):
    try:
        destination = generate_dashboard(args.input_dir, args.output)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(f"Dashboard saved to: {destination}")
    if not args.no_open:
        if not webbrowser.open(destination.as_uri()):
            print("No browser opened. Open the saved HTML file in your browser.")
