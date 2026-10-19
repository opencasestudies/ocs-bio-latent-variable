#!/usr/bin/env python3
"""Read-only, route-aware checks; never copy files, install software, or fit models."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROUTES = {
    "learner-r": ({"prepared", "config", "r_summary"}, set()),
    "learner-python": ({"prepared", "config", "python_summary"}, set()),
    "preprocess": ({"source"}, {"prepared", "config", "dense"}),
    "fit-r": ({"prepared"}, {"r_model", "checkpoint"}),
    "export-r": ({"prepared", "r_model"}, {"r_summary"}),
    "fit-python": ({"prepared", "dense"}, {"python_model"}),
    "inspect-python": ({"prepared", "python_model"}, set()),
}
GENERATED_INPUTS = {"prepared", "config", "dense", "r_model", "python_model"}
ROUTE_LABELS = {
    "learner-r": "Check the main learner files (R)",
    "learner-python": "Check the main learner files (Python)",
    "preprocess": "Rebuild preprocessing",
    "fit-r": "Fit the R model",
    "export-r": "Regenerate R exports",
    "fit-python": "Fit the Python model",
    "inspect-python": "Inspect the saved Python model",
}


def read_specs(project):
    manifest = json.loads((project / "data/external/reproduction_files.json").read_text())
    specs = {entry["id"]: entry for entry in manifest["files"]}
    for implementation in ("r", "python"):
        specs[f"{implementation}_summary"] = {
            "target_path": f"data/processed/selected_model_k6/{implementation}/pattern_summary.csv"
        }
    return specs


def inspect_file(path, expected=None):
    """Read all bytes; check reference identity only for original archive inputs."""
    try:
        if not path.exists():
            return False, "file is missing"
        if not path.is_file():
            return False, "path is not a regular file"
        size = path.stat().st_size
        if size == 0:
            return False, "file is empty"
        if expected and size != expected["size_bytes"]:
            return False, f"byte size mismatch: {size}; expected {expected['size_bytes']}"
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
                digest.update(block)
        if expected and digest.hexdigest() != expected["sha256"]:
            return False, "SHA-256 mismatch; do not use this as the archived input"
        suffix = "; SHA-256 matches archived reference" if expected else "; readable and nonempty"
        return True, f"{size:,} bytes{suffix}"
    except OSError as exc:
        return False, f"cannot read file: {exc}"


def check_route(project, archive, route, generated=(), specs=None):
    required, outputs = ROUTES[route]
    generated = set(generated)
    if not generated <= (required & GENERATED_INPUTS):
        raise ValueError("--generated-input must name a required, regenerable input for this route; source is never exempt.")
    specs = read_specs(project) if specs is None else specs
    rows = []
    for key, spec in specs.items():
        target = project / spec["target_path"]
        row = {"file": key, "path": str(target), "blocking": False}
        if key in outputs:
            row.update(status="output to be generated", detail="Produced by this stage; no download required.")
            if target.exists():
                row["detail"] = "Already exists. Work in a separate reproduction copy; this stage writes to this path."
                if key in {"r_model", "python_model"}:
                    row["blocking"] = True
                    row["detail"] += " The model launcher refuses to overwrite an existing full model."
        elif key not in required:
            row.update(status="not needed for this route", detail="Not a prerequisite for the selected stage.")
        else:
            expected = spec if "sha256" in spec and key not in generated else None
            ok, detail = inspect_file(target, expected)
            if ok:
                if key in generated:
                    detail += "; declared regenerated, not checked against original archive identity"
                row.update(status="ready", detail=detail)
            else:
                row.update(status="invalid required input" if target.exists() else "missing required input",
                           detail=detail, blocking=True)
                if "upload_filename" in spec and key not in generated:
                    archived = archive / spec["upload_filename"]
                    archive_ok, archive_detail = inspect_file(archived, spec)
                    row["archive_path"] = str(archived)
                    row["archive_check"] = archive_detail
                    if archive_ok:
                        row["detail"] += (
                            f"; verified archive copy available at {archived}. "
                            f"Copy it to {target}, then rerun this check. No files were copied."
                        )
        rows.append(row)

    try:
        if not archive.is_dir():
            archive_status = "archive directory missing"
        elif not any(path.is_file() for path in archive.iterdir()):
            archive_status = "archive directory empty"
        else:
            archive_status = "archive directory contains files"
    except OSError:
        archive_status = "archive directory unreadable"
    return {
        "route": route,
        "archive": {"path": str(archive), "status": archive_status,
                    "detail": "Folder existence does not download or verify inputs. Only files required for your route must be supplied."},
        "files": rows,
        "ready": not any(row["blocking"] for row in rows),
        "scope": "File prerequisites only; not environment, H5AD/RDS content, fitting, or scientific validation.",
    }


def format_summary(report, project):
    """Summarize the existing report without rechecking or changing any files."""
    def display_path(value):
        path = Path(value)
        try:
            return str(path.relative_to(project))
        except ValueError:
            return str(path)

    missing = [row for row in report["files"] if row["status"] == "missing required input"]
    invalid = [row for row in report["files"] if row["status"] == "invalid required input"]
    blocked_outputs = [row for row in report["files"]
                       if row["status"] == "output to be generated" and row["blocking"]]
    # The unchanged JSON report records existing outputs in its detail message.
    replacements = [row for row in report["files"]
                    if row["status"] == "output to be generated" and not row["blocking"]
                    and row["detail"].startswith("Already exists.")]
    lines = [f"Stage: {ROUTE_LABELS[report['route']]}"]
    if missing or invalid or blocked_outputs or replacements:
        lines.append("Paths below are relative to your project folder unless shown as absolute paths.")

    for heading, rows in (
        ("Missing required input:" if len(missing) == 1 else "Missing required inputs:", missing),
        ("Invalid required inputs:", invalid),
        ("Existing model files block fitting:", blocked_outputs),
        ("Files this stage will replace:", replacements),
    ):
        if not rows:
            continue
        lines.extend(["", heading])
        for row in rows:
            lines.append(f"  {display_path(row['path'])}")
            if row["status"] == "invalid required input":
                reason = row["detail"].split("; verified archive copy available at ", 1)[0]
                lines.append(f"    {reason}")
            if "archive_check" in row:
                archive_check = row["archive_check"]
                if "SHA-256 matches archived reference" in archive_check:
                    lines.append(f"    Verified download: {display_path(row['archive_path'])}")
                    lines.append("    Copy it to the project path above, then rerun this check.")
                elif archive_check != "file is missing":
                    lines.append(f"    Download check failed: {display_path(row['archive_path'])}")
                    lines.append(f"    {archive_check}")
        if rows is blocked_outputs:
            lines.append("  The model launcher refuses to overwrite these files. Use a separate working project.")

    if report["ready"]:
        conclusion = "Ready: the required files passed the checks."
    elif missing and not invalid and not blocked_outputs:
        noun = "input" if len(missing) == 1 else "inputs"
        conclusion = f"Not ready. Add the missing {noun} and run this check again."
    else:
        conclusion = "Not ready. Resolve the issues above and run this check again."
    lines.extend(["", conclusion, "This check has not changed any files."])
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--route", choices=ROUTES, required=True)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--archive", type=Path, help="Flat archive folder; defaults to the mounted project archive.")
    parser.add_argument("--generated-input", action="append", default=[], choices=sorted(GENERATED_INPUTS),
                        help="Explicitly identify an input you regenerated and checked in a preceding stage; may be repeated.")
    parser.add_argument("--json", action="store_true", help="Print the full machine-readable report instead of the brief summary.")
    args = parser.parse_args(argv)
    project = args.project_root.resolve()
    if not (project / "index.qmd").is_file() or not all((project / name).is_dir() for name in ("scripts", "data")):
        parser.error("Project root must contain index.qmd, scripts/, and data/.")
    archive = args.archive.resolve() if args.archive else project / "data/external/reproduction_archive"
    try:
        result = check_route(project, archive, args.route, args.generated_input)
    except (OSError, ValueError, KeyError) as exc:
        parser.error(str(exc))
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(format_summary(result, project))
    return 0 if result["ready"] else 1


if __name__ == "__main__":
    sys.exit(main())
