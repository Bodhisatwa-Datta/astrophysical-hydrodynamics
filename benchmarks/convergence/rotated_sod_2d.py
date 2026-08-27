"""Validate the 2D solver with uniform flow and rotated Sod tubes."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from time import perf_counter

import matplotlib
import numpy as np

from hydro.diagnostics import totals_2d
from hydro.exact_riemann import exact_riemann_solution
from hydro.solver2d import Grid2D, Solver2D

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def sod_x(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    left = x < 0.5
    return np.stack(
        (
            np.where(left, 1.0, 0.125),
            np.zeros_like(x),
            np.zeros_like(x),
            np.where(left, 1.0, 0.1),
        ),
        axis=0,
    )


def sod_y(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    lower = y < 0.5
    return np.stack(
        (
            np.where(lower, 1.0, 0.125),
            np.zeros_like(y),
            np.zeros_like(y),
            np.where(lower, 1.0, 0.1),
        ),
        axis=0,
    )


def uniform_validation(
    reconstruction: str, limiter: str, integrator: str
) -> dict[str, float | int]:
    grid = Grid2D(0.0, 1.0, 48, 0.0, 1.0, 36)
    solver = Solver2D(
        grid,
        cfl=0.4,
        riemann_solver="hllc",
        x_boundary="periodic",
        y_boundary="periodic",
        reconstruction=reconstruction,
        limiter=limiter,
        integrator=integrator,
    )
    primitive = np.empty((4, grid.nx, grid.ny))
    primitive[0] = 1.1
    primitive[1] = 0.3
    primitive[2] = -0.2
    primitive[3] = 0.9
    solver.initialise(primitive)
    initial = solver.active_conserved.copy()
    initial_totals = totals_2d(initial, grid.dx * grid.dy)
    solver.run(0.2)
    final_totals = totals_2d(solver.active_conserved, grid.dx * grid.dy)
    report: dict[str, float | int] = {
        "nx": grid.nx,
        "ny": grid.ny,
        "steps": solver.steps,
        "maximum_conserved_change": float(
            np.max(np.abs(solver.active_conserved - initial))
        ),
    }
    for quantity in initial_totals:
        report[f"relative_{quantity}_change"] = (
            final_totals[quantity] - initial_totals[quantity]
        ) / abs(initial_totals[quantity])
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[32, 64, 128])
    parser.add_argument("--time", type=float, default=0.1)
    parser.add_argument("--cfl", type=float, default=0.4)
    parser.add_argument(
        "--reconstruction", choices=("constant", "muscl"), default="constant"
    )
    parser.add_argument("--limiter", choices=("minmod", "mc", "vanleer"), default="mc")
    parser.add_argument("--integrator", choices=("euler", "rk2"), default="euler")
    parser.add_argument(
        "--csv", type=Path, default=Path("benchmarks/convergence/rotated_sod_2d.csv")
    )
    parser.add_argument(
        "--figure", type=Path, default=Path("figures/rotated_sod_2d_validation.png")
    )
    args = parser.parse_args()
    resolutions = sorted(set(args.resolutions))
    if len(resolutions) < 2 or resolutions[0] < 16:
        raise ValueError("provide at least two resolutions of 16 cells or more")

    rows: list[dict[str, float | int | str]] = []
    finest: tuple[Grid2D, np.ndarray, np.ndarray, np.ndarray] | None = None
    previous_error: float | None = None
    previous_resolution: int | None = None
    for resolution in resolutions:
        grid = Grid2D(0.0, 1.0, resolution, 0.0, 1.0, resolution)
        configuration = dict(
            cfl=args.cfl,
            riemann_solver="hll",
            reconstruction=args.reconstruction,
            limiter=args.limiter,
            integrator=args.integrator,
        )
        solver_x = Solver2D(grid, **configuration)
        solver_y = Solver2D(grid, **configuration)
        solver_x.initialise_function(sod_x)
        solver_y.initialise_function(sod_y)
        initial_totals = totals_2d(solver_x.active_conserved, grid.dx * grid.dy)
        start = perf_counter()
        solver_x.run(args.time)
        solver_y.run(args.time)
        runtime = perf_counter() - start
        primitive_x = solver_x.primitive
        primitive_y = solver_y.primitive
        exact = exact_riemann_solution(
            grid.x_centers,
            args.time,
            left=(1.0, 0.0, 1.0),
            right=(0.125, 0.0, 0.1),
            gamma=1.4,
            discontinuity=0.5,
        )
        density_x = primitive_x[0, :, grid.ny // 2]
        density_y = primitive_y[0, grid.nx // 2, :]
        error_x = float(grid.dx * np.sum(np.abs(density_x - exact[0])))
        error_y = float(grid.dy * np.sum(np.abs(density_y - exact[0])))
        order: float | str = ""
        if previous_error is not None and previous_resolution is not None:
            order = float(
                np.log(previous_error / error_x)
                / np.log(resolution / previous_resolution)
            )
        final_totals = totals_2d(solver_x.active_conserved, grid.dx * grid.dy)
        row: dict[str, float | int | str] = {
            "resolution": resolution,
            "reconstruction": args.reconstruction,
            "limiter": args.limiter if args.reconstruction == "muscl" else "none",
            "integrator": args.integrator,
            "steps": solver_x.steps,
            "runtime_two_runs_seconds": runtime,
            "density_L1_x": error_x,
            "density_L1_y": error_y,
            "observed_order": order,
            "rotational_density_Linf": float(
                np.max(np.abs(primitive_x[0] - primitive_y[0].T))
            ),
            "rotational_velocity_Linf": float(
                max(
                    np.max(np.abs(primitive_x[1] - primitive_y[2].T)),
                    np.max(np.abs(primitive_x[2] - primitive_y[1].T)),
                )
            ),
            "relative_mass_change": (
                final_totals["mass"] - initial_totals["mass"]
            )
            / initial_totals["mass"],
            "relative_energy_change": (
                final_totals["energy"] - initial_totals["energy"]
            )
            / initial_totals["energy"],
            "minimum_density": float(np.min(primitive_x[0])),
            "minimum_pressure": float(np.min(primitive_x[3])),
        }
        rows.append(row)
        previous_error = error_x
        previous_resolution = resolution
        finest = (grid, primitive_x, primitive_y, exact)

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    assert finest is not None
    grid, primitive_x, primitive_y, exact = finest
    figure, axes = plt.subplots(2, 2, figsize=(10.0, 7.5))
    image_x = axes[0, 0].pcolormesh(
        grid.x_centers,
        grid.y_centers,
        primitive_x[0].T,
        shading="auto",
        vmin=0.1,
        vmax=1.0,
    )
    axes[0, 0].set_title("Sod discontinuity normal to x")
    axes[0, 0].set_xlabel("x")
    axes[0, 0].set_ylabel("y")
    figure.colorbar(image_x, ax=axes[0, 0], label=r"Density $\rho$")
    image_y = axes[0, 1].pcolormesh(
        grid.x_centers,
        grid.y_centers,
        primitive_y[0].T,
        shading="auto",
        vmin=0.1,
        vmax=1.0,
    )
    axes[0, 1].set_title("Rotated discontinuity normal to y")
    axes[0, 1].set_xlabel("x")
    axes[0, 1].set_ylabel("y")
    figure.colorbar(image_y, ax=axes[0, 1], label=r"Density $\rho$")

    axes[1, 0].plot(grid.x_centers, exact[0], "k-", label="Exact")
    axes[1, 0].plot(
        grid.x_centers, primitive_x[0, :, grid.ny // 2], label="x-directed"
    )
    axes[1, 0].plot(
        grid.y_centers,
        primitive_y[0, grid.nx // 2, :],
        "--",
        label="y-directed",
    )
    axes[1, 0].set_xlabel("Normal coordinate")
    axes[1, 0].set_ylabel(r"Density $\rho$")
    axes[1, 0].set_title(f"Centre-line profiles, N={grid.nx}")
    axes[1, 0].grid(alpha=0.2)
    axes[1, 0].legend(frameon=False)

    axes[1, 1].loglog(
        resolutions,
        [float(row["density_L1_x"]) for row in rows],
        "o-",
        label="x-directed",
    )
    axes[1, 1].loglog(
        resolutions,
        [float(row["density_L1_y"]) for row in rows],
        "s--",
        label="y-directed",
    )
    axes[1, 1].set_xlabel("Cells per direction N")
    axes[1, 1].set_ylabel(r"Centre-line density $L_1$ error")
    axes[1, 1].set_title("Resolution dependence")
    axes[1, 1].grid(which="both", alpha=0.2)
    axes[1, 1].legend(frameon=False)
    method_name = (
        "first-order"
        if args.reconstruction == "constant" and args.integrator == "euler"
        else f"{args.reconstruction.upper()}-{args.limiter.upper()}/{args.integrator.upper()}"
    )
    figure.suptitle(f"{method_name} unsplit 2D HLL validation at t={args.time}")
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=180, bbox_inches="tight")

    print(
        json.dumps(
            {
                "uniform_flow": uniform_validation(
                    args.reconstruction, args.limiter, args.integrator
                ),
                "rotated_sod": rows,
            },
            indent=2,
        )
    )
    print(f"Saved {args.csv}")
    print(f"Saved {args.figure}")


if __name__ == "__main__":
    main()

