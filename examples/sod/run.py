"""Run the first-order HLL solver on the Sod shock tube."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from hydro.diagnostics import state_summary, totals
from hydro.exact_riemann import exact_riemann_solution
from hydro.problems import sod_initial_condition
from hydro.solver1d import Grid1D, Solver1D


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, default=400)
    parser.add_argument("--time", type=float, default=0.2)
    parser.add_argument("--cfl", type=float, default=0.8)
    parser.add_argument("--output", type=Path, default=Path("figures/sod_hll_first_order.png"))
    args = parser.parse_args()

    gamma = 1.4
    grid = Grid1D(0.0, 1.0, args.cells)
    solver = Solver1D(grid, gamma=gamma, cfl=args.cfl)
    solver.initialise_function(sod_initial_condition)
    initial_totals = totals(solver.active_conserved, grid.dx)
    solver.run(args.time)
    numerical = solver.primitive
    exact = exact_riemann_solution(
        grid.centers,
        args.time,
        left=(1.0, 0.0, 1.0),
        right=(0.125, 0.0, 0.1),
        gamma=gamma,
        discontinuity=0.5,
    )

    numerical_internal = numerical[2] / ((gamma - 1.0) * numerical[0])
    exact_internal = exact[2] / ((gamma - 1.0) * exact[0])
    fields = (
        (numerical[0], exact[0], r"Density $\rho$"),
        (numerical[1], exact[1], r"Velocity $u$"),
        (numerical[2], exact[2], r"Pressure $p$"),
        (numerical_internal, exact_internal, r"Specific internal energy $e$"),
    )
    figure, axes = plt.subplots(2, 2, figsize=(9.0, 6.5), sharex=True)
    for axis, (values, reference, ylabel) in zip(axes.flat, fields):
        axis.plot(grid.centers, reference, "k-", linewidth=1.3, label="Exact")
        axis.plot(grid.centers, values, color="tab:blue", linewidth=1.0, label="HLL")
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.2)
    axes[0, 0].legend(frameon=False)
    axes[1, 0].set_xlabel("Position $x$ (dimensionless)")
    axes[1, 1].set_xlabel("Position $x$ (dimensionless)")
    figure.suptitle(f"Sod shock tube, t={solver.time:.3f}, N={grid.n_cells}")
    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180, bbox_inches="tight")

    final_totals = totals(solver.active_conserved, grid.dx)
    report = {
        "cells": grid.n_cells,
        "time": solver.time,
        "steps": solver.steps,
        "density_L1_error": float(grid.dx * np.sum(np.abs(numerical[0] - exact[0]))),
        "relative_mass_change": (final_totals["mass"] - initial_totals["mass"]) / initial_totals["mass"],
        "relative_energy_change": (final_totals["energy"] - initial_totals["energy"]) / initial_totals["energy"],
        **state_summary(solver.active_conserved, gamma),
    }
    print(json.dumps(report, indent=2))
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()

