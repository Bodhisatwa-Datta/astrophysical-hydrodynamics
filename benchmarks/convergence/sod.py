"""Measure resolution-dependent errors for the first-order Sod solution."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from hydro.problems import SOD
from hydro.validation import run_riemann_problem


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[50, 100, 200, 400, 800])
    parser.add_argument("--cfl", type=float, default=0.8)
    parser.add_argument(
        "--csv", type=Path, default=Path("benchmarks/convergence/sod_first_order.csv")
    )
    parser.add_argument(
        "--figure", type=Path, default=Path("figures/sod_convergence_first_order.png")
    )
    args = parser.parse_args()
    resolutions = sorted(set(args.resolutions))
    if len(resolutions) < 2 or resolutions[0] < 10:
        raise ValueError("provide at least two resolutions of 10 cells or more")

    rows: list[dict[str, float | int | str]] = []
    previous_errors: dict[str, float] = {}
    for cells in resolutions:
        result = run_riemann_problem(SOD, cells, args.cfl)
        row: dict[str, float | int | str] = {
            "cells": cells,
            "dx": 1.0 / cells,
            "steps": result.diagnostics["steps"],
            "runtime_seconds": result.diagnostics["runtime_seconds"],
        }
        for field in ("density", "velocity", "pressure"):
            error = float(result.diagnostics[f"{field}_L1_error"])
            row[f"{field}_L1_error"] = error
            if rows:
                refinement = cells / int(rows[-1]["cells"])
                row[f"{field}_observed_order"] = float(
                    np.log(previous_errors[field] / error) / np.log(refinement)
                )
            else:
                row[f"{field}_observed_order"] = ""
            previous_errors[field] = error
        rows.append(row)

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    figure, axis = plt.subplots(figsize=(6.5, 4.5))
    for field, marker in zip(("density", "velocity", "pressure"), ("o", "s", "^")):
        errors = [float(row[f"{field}_L1_error"]) for row in rows]
        axis.loglog(resolutions, errors, marker=marker, label=field.title())
    reference = float(rows[0]["density_L1_error"]) * resolutions[0] / np.asarray(resolutions)
    axis.loglog(resolutions, reference, "k--", alpha=0.7, label=r"$N^{-1}$ reference")
    axis.set_xlabel("Number of cells $N$")
    axis.set_ylabel(r"$L_1$ error")
    axis.set_title("Sod shock tube: first-order HLL resolution study")
    axis.grid(which="both", alpha=0.25)
    axis.legend(frameon=False)
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=180, bbox_inches="tight")

    print(json.dumps(rows, indent=2))
    print(f"Saved {args.csv}")
    print(f"Saved {args.figure}")


if __name__ == "__main__":
    main()
