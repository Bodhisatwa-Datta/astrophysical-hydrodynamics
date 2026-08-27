"""Run and plot any configured one-dimensional Riemann benchmark."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from hydro.problems import RIEMANN_PROBLEMS
from hydro.validation import contact_transition_width, run_riemann_problem


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("problem", choices=RIEMANN_PROBLEMS)
    parser.add_argument("--cells", type=int, default=400)
    parser.add_argument("--cfl", type=float, default=0.8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    problem = RIEMANN_PROBLEMS[args.problem]
    output = args.output or Path("figures") / f"{problem.name}_hll_first_order.png"
    result = run_riemann_problem(problem, args.cells, args.cfl)

    gamma = problem.gamma
    numerical_internal = result.numerical[2] / ((gamma - 1.0) * result.numerical[0])
    exact_internal = result.exact[2] / ((gamma - 1.0) * result.exact[0])
    fields = (
        (result.numerical[0], result.exact[0], r"Density $\rho$"),
        (result.numerical[1], result.exact[1], r"Velocity $u$"),
        (result.numerical[2], result.exact[2], r"Pressure $p$"),
        (numerical_internal, exact_internal, r"Specific internal energy $e$"),
    )
    figure, axes = plt.subplots(2, 2, figsize=(9.0, 6.5), sharex=True)
    for axis, (values, reference, ylabel) in zip(axes.flat, fields):
        axis.plot(result.x, reference, "k-", linewidth=1.3, label="Exact")
        axis.plot(result.x, values, color="tab:blue", linewidth=1.0, label="HLL")
        combined = np.concatenate((np.ravel(values), np.ravel(reference)))
        scale = max(float(np.max(np.abs(combined))), 1.0)
        if float(np.ptp(combined)) < 1.0e-10 * scale:
            centre = float(np.mean(combined))
            axis.set_ylim(centre - 0.05 * scale, centre + 0.05 * scale)
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.2)
    axes[0, 0].legend(frameon=False)
    axes[1, 0].set_xlabel("Position $x$ (dimensionless)")
    axes[1, 1].set_xlabel("Position $x$ (dimensionless)")
    display_name = problem.name.replace("_", " ").title()
    figure.suptitle(
        f"{display_name}, t={problem.final_time:g}, N={args.cells}"
    )
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180, bbox_inches="tight")

    report = dict(result.diagnostics)
    if problem.name == "contact_discontinuity":
        report["contact_transition_width"] = contact_transition_width(result)
    print(json.dumps(report, indent=2))
    print(f"Saved {output}")


if __name__ == "__main__":
    main()
