"""Compare first-order and limited MUSCL solutions on four Riemann problems."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt

from hydro.problems import CONTACT_DISCONTINUITY, SOD, STRONG_RAREFACTION, STRONG_SHOCK
from hydro.validation import contact_transition_width, run_riemann_problem


METHODS = (
    ("first_order", "constant", "mc", "euler"),
    ("muscl_minmod", "muscl", "minmod", "rk2"),
    ("muscl_mc", "muscl", "mc", "rk2"),
    ("muscl_vanleer", "muscl", "vanleer", "rk2"),
)
PROBLEMS = (SOD, CONTACT_DISCONTINUITY, STRONG_SHOCK, STRONG_RAREFACTION)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, default=400)
    parser.add_argument("--cfl", type=float, default=0.4)
    parser.add_argument(
        "--csv", type=Path, default=Path("benchmarks/reconstruction/comparison.csv")
    )
    parser.add_argument(
        "--figure", type=Path, default=Path("figures/reconstruction_comparison.png")
    )
    args = parser.parse_args()

    rows: list[dict[str, float | int | str]] = []
    results = {}
    for problem in PROBLEMS:
        for label, reconstruction, limiter, integrator in METHODS:
            result = run_riemann_problem(
                problem,
                args.cells,
                args.cfl,
                reconstruction,
                limiter,
                integrator,
            )
            results[(problem.name, label)] = result
            row: dict[str, float | int | str] = {
                "problem": problem.name,
                "method": label,
                "cells": args.cells,
                "steps": result.diagnostics["steps"],
                "runtime_seconds": result.diagnostics["runtime_seconds"],
                "density_L1_error": result.diagnostics["density_L1_error"],
                "pressure_L1_error": result.diagnostics["pressure_L1_error"],
                "minimum_density": result.diagnostics["minimum_density"],
                "minimum_pressure": result.diagnostics["minimum_pressure"],
                "contact_transition_width": "",
            }
            if problem.name == "contact_discontinuity":
                row["contact_transition_width"] = contact_transition_width(result)
            rows.append(row)

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    figure, axes = plt.subplots(2, 2, figsize=(10.0, 7.2), sharex=True)
    colors = ("tab:gray", "tab:orange", "tab:blue", "tab:green")
    for axis, problem in zip(axes.flat, PROBLEMS):
        first = results[(problem.name, METHODS[0][0])]
        axis.plot(first.x, first.exact[0], "k-", linewidth=1.5, label="Exact")
        for (label, *_), color in zip(METHODS, colors):
            result = results[(problem.name, label)]
            axis.plot(
                result.x,
                result.numerical[0],
                color=color,
                linewidth=1.0,
                label=label.replace("_", " ").title(),
            )
        axis.set_title(problem.name.replace("_", " ").title())
        axis.set_ylabel(r"Density $\rho$")
        axis.grid(alpha=0.2)
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[1, 0].set_xlabel("Position $x$ (dimensionless)")
    axes[1, 1].set_xlabel("Position $x$ (dimensionless)")
    figure.suptitle(f"Reconstruction comparison, N={args.cells}, CFL={args.cfl}")
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=180, bbox_inches="tight")

    print(json.dumps(rows, indent=2))
    print(f"Saved {args.csv}")
    print(f"Saved {args.figure}")


if __name__ == "__main__":
    main()

