"""Measure smooth-wave accuracy of first-order and MUSCL-RK2 methods."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from hydro.validation import run_entropy_wave


METHODS = (
    ("first_order", "constant", "mc", "euler"),
    ("muscl_minmod", "muscl", "minmod", "rk2"),
    ("muscl_mc", "muscl", "mc", "rk2"),
    ("muscl_vanleer", "muscl", "vanleer", "rk2"),
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[32, 64, 128, 256])
    parser.add_argument("--cfl", type=float, default=0.4)
    parser.add_argument(
        "--csv", type=Path, default=Path("benchmarks/convergence/entropy_wave.csv")
    )
    parser.add_argument(
        "--figure", type=Path, default=Path("figures/entropy_wave_convergence.png")
    )
    args = parser.parse_args()
    resolutions = sorted(set(args.resolutions))
    if len(resolutions) < 2 or resolutions[0] < 16:
        raise ValueError("provide at least two resolutions of 16 cells or more")

    rows: list[dict[str, float | int | str]] = []
    for label, reconstruction, limiter, integrator in METHODS:
        previous_cells: int | None = None
        previous_error: float | None = None
        for cells in resolutions:
            result = run_entropy_wave(
                cells, args.cfl, reconstruction, limiter, integrator
            )
            error = float(result.diagnostics["density_L1_error"])
            order: float | str = ""
            if previous_error is not None and previous_cells is not None:
                order = float(
                    np.log(previous_error / error)
                    / np.log(cells / previous_cells)
                )
            rows.append(
                {
                    "method": label,
                    "cells": cells,
                    "dx": 1.0 / cells,
                    "steps": result.diagnostics["steps"],
                    "runtime_seconds": result.diagnostics["runtime_seconds"],
                    "density_L1_error": error,
                    "observed_order": order,
                    "relative_mass_change": result.diagnostics["relative_mass_change"],
                }
            )
            previous_cells = cells
            previous_error = error

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    figure, axis = plt.subplots(figsize=(7.0, 4.8))
    markers = ("o", "s", "^", "D")
    for (label, *_), marker in zip(METHODS, markers):
        selected = [row for row in rows if row["method"] == label]
        axis.loglog(
            [int(row["cells"]) for row in selected],
            [float(row["density_L1_error"]) for row in selected],
            marker=marker,
            label=label.replace("_", " ").title(),
        )
    reference_n = np.asarray(resolutions, dtype=float)
    finest_mc = next(
        float(row["density_L1_error"])
        for row in reversed(rows)
        if row["method"] == "muscl_mc"
    )
    second_order_reference = finest_mc * (resolutions[-1] / reference_n) ** 2
    axis.loglog(reference_n, second_order_reference, "k--", alpha=0.6, label=r"$N^{-2}$")
    axis.set_xlabel("Number of cells $N$")
    axis.set_ylabel(r"Density $L_1$ error")
    axis.set_title("Periodic entropy wave after one crossing")
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

