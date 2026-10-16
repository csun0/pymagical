"""Compare supplied candidate peaks with unique peaks in reported circuits.

Run from the repository with: uv run eval/analysis/plot_peak_input_output.py
"""
from plot_gene_input_output import plot_input_output


if __name__ == "__main__":
    plot_input_output("peaks")
