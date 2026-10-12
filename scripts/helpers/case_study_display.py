"""Keep learner output readable in both a Terminal and the Quarto page."""

import os
from pathlib import Path

import matplotlib.pyplot as plt

from .case_study_paths import ROOT


def _in_quarto_render():
    # Quarto sets this for executed cells; a plain Python session does not.
    return bool(os.environ.get("QUARTO_DOCUMENT_PATH"))


def render_table(data_frame, *, classes="table table-sm table-striped", escape=False):
    """Return HTML during rendering, otherwise an untruncated plain-text table."""
    if _in_quarto_render():
        return data_frame.to_html(
            index=False, escape=escape, border=0, classes=classes
        )
    return data_frame.to_string(index=False, line_width=100)


def show_plot(figure, filename):
    """Display a plot in Quarto or save a learner PNG without changing the backend."""
    if _in_quarto_render():
        plt.show()
        return

    if Path(filename).name != filename or Path(filename).suffix != ".png":
        raise ValueError("Use a PNG filename without a directory path.")
    output_dir = ROOT / "learner_outputs" / "python"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / filename
    figure.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(figure)
    print(f"Saved plot: {output_path}")
