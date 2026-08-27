"""Compare HLL, HLLC, and Rusanov fluxes under matched discretizations."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt

from hydro.problems import CONTACT_DISCONTINUITY, SOD, STRONG_RAREFACTION, STRONG_SHOCK
from hydro.validation import contact_transition_width, run_riemann_problem


FLUXES = ("hll", "hllc", "rusanov")
CONFIGURATIONS = (
    ("first_order", "constant", "mc", "euler"),
    ("muscl_mc", "muscl", "mc", "rk2"),
)
PROBLEMS = (SOD, CONTACT_DISCONTINUITY, STRONG_SHOCK, STRONG_RAREFACTION)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, default=400)
    parser.add_argument("--cfl", type=float, default=0.4)
    parser.add_argument(
        "--csv", type=Path, default=Path("benchmarks/riemann_solvers/comparison.csv")
    )
    parser.add_argument(
        "--figure", type=Path, default=Path("figures/riemann_solver_comparison.png")
    )
    args = parser.parse_args()

    rows: list[dict[str, float | int | str]] = []
    results = {}
    for problem in PROBLEMS:
        for configuration, reconstruction, limiter, integrator in CONFIGURATIONS:
            for flux in FLUXES:
                result = run_riemann_problem(
                    problem,
                    args.cells,
                    args.cfl,
                    reconstruction,
                    limiter,
                    integrator,
                    flux,
                )
                results[(problem.name, configuration, flux)] = result
                row: dict[str, float | int | str] = {
                    "problem": problem.name,
                    "configuration": configuration,
                    "flux": flux,
                    "cells": args.cells,
                    "steps": result.diagnostics["steps"],
                    "runtime_seconds": result.diagnostics["runtime_seconds"],
                    "density_L1_error": result.diagnostics["density_L1_error"],
                    "pressure_L1_error": result.diagnostics["pressure_L1_error"],
                    "minimum_density": result.diagnostics["minimum_density"],
                    "minimum_pressure": result.diagnostics["minimum_pressure"],
                    "relative_mass_budget_residual": result.diagnostics[
                        "relative_mass_boundary_budget_residual"
                    ],
                    "relative_energy_budget_residual": result.diagnostics[
                        "relative_energy_boundary_budget_residual"
                    ],
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

    figure, axes = plt.subplots(2, 4, figsize=(14.0, 6.4), sharex="col")
    colors = {"hll": "tab:blue", "hllc": "tab:green", "rusanov": "tab:red"}
    for row_index, (configuration, *_methods) in enumerate(CONFIGURATIONS):
        for column_index, problem in enumerate(PROBLEMS):
            axis = axes[row_index, column_index]
            first = results[(problem.name, configuration, FLUXES[0])]
            axis.plot(first.x, first.exact[0], "k-", linewidth=1.4, label="Exact")
            for flux in FLUXES:
                result = results[(problem.name, configuration, flux)]
                axis.plot(
                    result.x,
                    result.numerical[0],
                    color=colors[flux],
                    linewidth=1.0,
                    label=flux.upper(),
                )
            if row_index == 0:
                axis.set_title(problem.name.replace("_", " ").title())
            if column_index == 0:
                axis.set_ylabel(
                    f"{configuration.replace('_', ' ').title()}\nDensity $\\rho$"
                )
            if row_index == 1:
                axis.set_xlabel("Position $x$")
            axis.grid(alpha=0.2)
    axes[0, 0].legend(frameon=False, fontsize=8)
    figure.suptitle(
        f"Riemann-flux comparison, N={args.cells}, CFL={args.cfl}"
    )
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=180, bbox_inches="tight")

    print(json.dumps(rows, indent=2))
    print(f"Saved {args.csv}")
    print(f"Saved {args.figure}")


if __name__ == "__main__":
    main()

