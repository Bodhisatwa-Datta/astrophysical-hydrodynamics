"""Measure convergence for diagonal entropy-wave advection in two dimensions."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib
import numpy as np

from hydro.validation import run_entropy_wave_2d

matplotlib.use("Agg")
import matplotlib.pyplot as plt


METHODS = {
    "first_order": ("constant", "mc", "euler"),
    "muscl_minmod": ("muscl", "minmod", "rk2"),
    "muscl_mc": ("muscl", "mc", "rk2"),
    "muscl_vanleer": ("muscl", "vanleer", "rk2"),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[16, 32, 64])
    parser.add_argument("--cfl", type=float, default=0.4)
    parser.add_argument(
        "--csv", type=Path, default=Path("benchmarks/convergence/entropy_wave_2d.csv")
    )
    parser.add_argument(
        "--figure", type=Path, default=Path("figures/entropy_wave_2d_convergence.png")
    )
    args = parser.parse_args()
    resolutions = sorted(set(args.resolutions))
    if len(resolutions) < 2 or resolutions[0] < 8:
        raise ValueError("provide at least two resolutions of 8 cells or more")

    rows: list[dict[str, float | int | str]] = []
    for method, (reconstruction, limiter, integrator) in METHODS.items():
        previous_error: float | None = None
        previous_resolution: int | None = None
        for resolution in resolutions:
            result = run_entropy_wave_2d(
                resolution,
                cfl=args.cfl,
                reconstruction=reconstruction,
                limiter=limiter,
                integrator=integrator,
            )
            error = float(result.diagnostics["density_L1_error"])
            order: float | str = ""
            if previous_error is not None and previous_resolution is not None:
                order = float(
                    np.log(previous_error / error)
                    / np.log(resolution / previous_resolution)
                )
            rows.append({"method": method, "observed_order": order, **result.diagnostics})
            previous_error = error
            previous_resolution = resolution

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    figure, axis = plt.subplots(figsize=(7.2, 5.2))
    for method in METHODS:
        selected = [row for row in rows if row["method"] == method]
        axis.loglog(
            [int(row["resolution"]) for row in selected],
            [float(row["density_L1_error"]) for row in selected],
            "o-",
            label=method.replace("_", " "),
        )
    reference_resolution = np.asarray(resolutions, dtype=float)
    finest_mc = next(
        float(row["density_L1_error"])
        for row in rows
        if row["method"] == "muscl_mc" and row["resolution"] == resolutions[-1]
    )
    reference = finest_mc * (resolutions[-1] / reference_resolution) ** 2
    axis.loglog(reference_resolution, reference, "k:", label=r"$N^{-2}$ reference")
    axis.set_xlabel("Cells per direction N")
    axis.set_ylabel(r"Density $L_1$ error")
    axis.set_title("Diagonal periodic entropy-wave convergence")
    axis.grid(which="both", alpha=0.2)
    axis.legend(frameon=False)
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=180, bbox_inches="tight")

    print(json.dumps(rows, indent=2))
    print(f"Saved {args.csv}")
    print(f"Saved {args.figure}")


if __name__ == "__main__":
    main()
