"""Run a smooth periodic Kelvin--Helmholtz resolution study."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from time import perf_counter

import matplotlib
import numpy as np

from hydro.diagnostics import (
    state_summary_2d,
    totals_2d,
    transverse_kinetic_energy_2d,
    vorticity_2d,
)
from hydro.problems import kelvin_helmholtz_2d
from hydro.solver2d import Grid2D, Solver2D

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def record_diagnostics(
    solver: Solver2D,
    initial_totals: dict[str, float],
    runtime: float,
) -> dict[str, float | int]:
    grid = solver.grid
    primitive = solver.primitive
    cell_area = grid.dx * grid.dy
    totals = totals_2d(solver.active_conserved, cell_area)
    vorticity = vorticity_2d(primitive, grid.dx, grid.dy)
    record: dict[str, float | int] = {
        "resolution": grid.nx,
        "time": solver.time,
        "steps": solver.steps,
        "runtime_seconds": runtime,
        "transverse_kinetic_energy": transverse_kinetic_energy_2d(
            primitive, cell_area
        ),
        "rms_vertical_velocity": float(np.sqrt(np.mean(primitive[2] ** 2))),
        "rms_vorticity": float(np.sqrt(np.mean(vorticity**2))),
        "maximum_absolute_vorticity": float(np.max(np.abs(vorticity))),
        "density_standard_deviation": float(np.std(primitive[0])),
    }
    for quantity, initial in initial_totals.items():
        change = totals[quantity] - initial
        record[f"absolute_{quantity}_change"] = change
        if quantity != "momentum_y":
            record[f"relative_{quantity}_change"] = change / abs(initial)
    record.update(state_summary_2d(solver.active_conserved, solver.gamma))
    return record


def run_resolution(
    resolution: int,
    final_time: float,
    sample_interval: float,
    cfl: float,
) -> tuple[list[dict[str, float | int]], np.ndarray, Grid2D]:
    grid = Grid2D(0.0, 1.0, resolution, 0.0, 1.0, resolution)
    solver = Solver2D(
        grid,
        gamma=1.4,
        cfl=cfl,
        riemann_solver="hllc",
        reconstruction="muscl",
        limiter="mc",
        integrator="rk2",
        x_boundary="periodic",
        y_boundary="periodic",
    )
    solver.initialise_function(kelvin_helmholtz_2d)
    initial_totals = totals_2d(solver.active_conserved, grid.dx * grid.dy)
    records = [record_diagnostics(solver, initial_totals, 0.0)]
    start = perf_counter()
    target = min(sample_interval, final_time)
    while solver.time < final_time:
        solver.step(min(solver.timestep(), target - solver.time))
        if np.isclose(solver.time, target, rtol=0.0, atol=32.0 * np.finfo(float).eps):
            records.append(record_diagnostics(solver, initial_totals, perf_counter() - start))
            target = min(target + sample_interval, final_time)
            if records[-1]["time"] == final_time:
                break
    return records, solver.primitive.copy(), grid


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[64, 96, 128])
    parser.add_argument("--time", type=float, default=1.5)
    parser.add_argument("--sample-interval", type=float, default=0.05)
    parser.add_argument("--cfl", type=float, default=0.35)
    parser.add_argument(
        "--history", type=Path, default=Path("benchmarks/kelvin_helmholtz_history.csv")
    )
    parser.add_argument(
        "--summary", type=Path, default=Path("benchmarks/kelvin_helmholtz_summary.csv")
    )
    parser.add_argument(
        "--figure", type=Path, default=Path("figures/kelvin_helmholtz_fields.png")
    )
    parser.add_argument(
        "--growth-figure",
        type=Path,
        default=Path("figures/kelvin_helmholtz_growth.png"),
    )
    args = parser.parse_args()
    resolutions = sorted(set(args.resolutions))
    if len(resolutions) < 2 or resolutions[0] < 32:
        raise ValueError("provide at least two resolutions of 32 cells or more")
    if args.time <= 0.0 or args.sample_interval <= 0.0:
        raise ValueError("time and sample interval must be positive")

    all_records: list[dict[str, float | int]] = []
    final_states: list[tuple[int, np.ndarray, Grid2D]] = []
    summaries: list[dict[str, float | int]] = []
    for resolution in resolutions:
        records, primitive, grid = run_resolution(
            resolution, args.time, args.sample_interval, args.cfl
        )
        all_records.extend(records)
        final_states.append((resolution, primitive, grid))
        initial = records[0]
        final = records[-1]
        summaries.append(
            {
                **final,
                "transverse_energy_growth": float(final["transverse_kinetic_energy"])
                / float(initial["transverse_kinetic_energy"]),
                "vertical_velocity_growth": float(final["rms_vertical_velocity"])
                / float(initial["rms_vertical_velocity"]),
            }
        )

    args.history.parent.mkdir(parents=True, exist_ok=True)
    with args.history.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(all_records[0]))
        writer.writeheader()
        writer.writerows(all_records)
    with args.summary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)

    figure, axes = plt.subplots(
        len(final_states), 3, figsize=(12.0, 3.6 * len(final_states)), squeeze=False
    )
    for row, (resolution, primitive, grid) in enumerate(final_states):
        vorticity = vorticity_2d(primitive, grid.dx, grid.dy)
        speed = np.sqrt(primitive[1] ** 2 + primitive[2] ** 2)
        fields = (primitive[0], vorticity, speed)
        labels = (r"Density $\rho$", r"Vorticity $\omega_z$", r"Speed $|\mathbf{v}|$")
        color_maps = ("viridis", "RdBu_r", "magma")
        for column, (field, label, color_map) in enumerate(
            zip(fields, labels, color_maps)
        ):
            limit = float(np.max(np.abs(field))) if column == 1 else None
            image = axes[row, column].pcolormesh(
                grid.x_centers,
                grid.y_centers,
                field.T,
                shading="auto",
                cmap=color_map,
                vmin=-limit if limit is not None else None,
                vmax=limit,
            )
            axes[row, column].set_aspect("equal")
            axes[row, column].set_xlabel("x")
            axes[row, column].set_ylabel("y")
            axes[row, column].set_title(f"{label}, {resolution}²")
            figure.colorbar(image, ax=axes[row, column], label=label)
    figure.suptitle(f"Smooth Kelvin--Helmholtz flow at t={args.time:g} (dimensionless)")
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=180, bbox_inches="tight")

    growth_figure, growth_axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    for resolution in resolutions:
        selected = [row for row in all_records if row["resolution"] == resolution]
        times = [float(row["time"]) for row in selected]
        initial_energy = float(selected[0]["transverse_kinetic_energy"])
        growth_axes[0].semilogy(
            times,
            [float(row["transverse_kinetic_energy"]) / initial_energy for row in selected],
            label=f"{resolution}²",
        )
        growth_axes[1].plot(
            times,
            [float(row["maximum_absolute_vorticity"]) for row in selected],
            label=f"{resolution}²",
        )
    growth_axes[0].set_ylabel("Transverse kinetic energy / initial value")
    growth_axes[1].set_ylabel(r"Maximum $|\omega_z|$")
    for axis in growth_axes:
        axis.set_xlabel("Time")
        axis.grid(alpha=0.2)
        axis.legend(frameon=False)
    growth_figure.suptitle("Kelvin--Helmholtz growth diagnostics")
    growth_figure.tight_layout()
    args.growth_figure.parent.mkdir(parents=True, exist_ok=True)
    growth_figure.savefig(args.growth_figure, dpi=180, bbox_inches="tight")

    print(json.dumps(summaries, indent=2))
    print(f"Saved {args.history}")
    print(f"Saved {args.summary}")
    print(f"Saved {args.figure}")
    print(f"Saved {args.growth_figure}")


if __name__ == "__main__":
    main()
