"""Run perturbed and hydrostatic-control Rayleigh--Taylor calculations."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from time import perf_counter

import matplotlib
import numpy as np

from hydro.diagnostics import (
    density_interface_height_2d,
    state_summary_2d,
    total_energy_with_gravity_2d,
    totals_2d,
    transverse_kinetic_energy_2d,
    vorticity_2d,
)
from hydro.problems import rayleigh_taylor_2d
from hydro.solver2d import Grid2D, Solver2D

matplotlib.use("Agg")
import matplotlib.pyplot as plt


GRAVITY = (0.0, -0.5)


def record(
    solver: Solver2D,
    case: str,
    initial_totals: dict[str, float],
    initial_total_energy: float,
    runtime: float,
) -> dict[str, float | int | str]:
    grid = solver.grid
    primitive = solver.primitive
    cell_area = grid.dx * grid.dy
    totals = totals_2d(solver.active_conserved, cell_area)
    total_energy = total_energy_with_gravity_2d(
        solver.active_conserved, *grid.mesh, solver.gravity, cell_area
    )
    heights = density_interface_height_2d(primitive, grid.y_centers)
    if np.any(~np.isfinite(heights)):
        raise RuntimeError("density interface could not be located in every column")
    diagnostics: dict[str, float | int | str] = {
        "case": case,
        "resolution": grid.nx,
        "time": solver.time,
        "steps": solver.steps,
        "runtime_seconds": runtime,
        "transverse_kinetic_energy": transverse_kinetic_energy_2d(
            primitive, cell_area
        ),
        "rms_vertical_velocity": float(np.sqrt(np.mean(primitive[2] ** 2))),
        "maximum_vertical_velocity": float(np.max(np.abs(primitive[2]))),
        "bubble_height": float(np.max(heights) - 0.5),
        "spike_depth": float(0.5 - np.min(heights)),
        "maximum_absolute_vorticity": float(
            np.max(
                np.abs(
                    vorticity_2d(
                        primitive, grid.dx, grid.dy, periodic=(True, False)
                    )
                )
            )
        ),
        "relative_gas_plus_potential_energy_change": (
            total_energy - initial_total_energy
        )
        / abs(initial_total_energy),
    }
    diagnostics["relative_mass_change"] = (
        totals["mass"] - initial_totals["mass"]
    ) / initial_totals["mass"]
    diagnostics["absolute_momentum_x_change"] = (
        totals["momentum_x"] - initial_totals["momentum_x"]
    )
    diagnostics["absolute_momentum_y_change"] = (
        totals["momentum_y"] - initial_totals["momentum_y"]
    )
    diagnostics.update(state_summary_2d(solver.active_conserved, solver.gamma))
    return diagnostics


def run_case(
    resolution: int,
    case: str,
    final_time: float,
    sample_interval: float,
    cfl: float,
) -> tuple[list[dict[str, float | int | str]], np.ndarray, Grid2D]:
    perturbation = 0.0025 if case == "perturbed" else 0.0
    grid = Grid2D(0.0, 1.0, resolution, 0.0, 1.0, resolution)
    solver = Solver2D(
        grid,
        gamma=1.4,
        cfl=cfl,
        riemann_solver="hllc",
        reconstruction="muscl",
        limiter="mc",
        integrator="rk2",
        gravity=GRAVITY,
        x_boundary="periodic",
        y_boundary="hydrostatic",
    )
    solver.initialise_function(
        lambda x, y: rayleigh_taylor_2d(
            x, y, gravity_y=GRAVITY[1], perturbation_amplitude=perturbation
        )
    )
    cell_area = grid.dx * grid.dy
    initial_totals = totals_2d(solver.active_conserved, cell_area)
    initial_total_energy = total_energy_with_gravity_2d(
        solver.active_conserved, *grid.mesh, GRAVITY, cell_area
    )
    records = [record(solver, case, initial_totals, initial_total_energy, 0.0)]
    start = perf_counter()
    target = min(sample_interval, final_time)
    while solver.time < final_time:
        solver.step(min(solver.timestep(), target - solver.time))
        if np.isclose(solver.time, target, rtol=0.0, atol=32.0 * np.finfo(float).eps):
            records.append(
                record(
                    solver,
                    case,
                    initial_totals,
                    initial_total_energy,
                    perf_counter() - start,
                )
            )
            target = min(target + sample_interval, final_time)
            if records[-1]["time"] == final_time:
                break
    return records, solver.primitive.copy(), grid


def measured_interface_growth_rate(
    records: list[dict[str, float | int | str]],
    start_time: float = 0.8,
    end_time: float = 2.2,
) -> float:
    """Fit exponential interface-amplitude growth over a stated time window."""
    end_time = min(end_time, float(records[-1]["time"]))
    start_time = min(start_time, 0.4 * end_time)
    selected = [
        row for row in records if start_time <= float(row["time"]) <= end_time
    ]
    times = np.asarray([float(row["time"]) for row in selected])
    amplitudes = np.asarray(
        [
            0.5 * (float(row["bubble_height"]) + float(row["spike_depth"]))
            for row in selected
        ]
    )
    if times.size < 2 or np.any(amplitudes <= 0.0):
        raise ValueError("growth-rate window requires two positive interface amplitudes")
    return float(np.polyfit(times, np.log(amplitudes), 1)[0])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[64, 96, 128])
    parser.add_argument("--time", type=float, default=2.5)
    parser.add_argument("--sample-interval", type=float, default=0.1)
    parser.add_argument("--cfl", type=float, default=0.3)
    parser.add_argument(
        "--history", type=Path, default=Path("benchmarks/rayleigh_taylor_history.csv")
    )
    parser.add_argument(
        "--summary", type=Path, default=Path("benchmarks/rayleigh_taylor_summary.csv")
    )
    parser.add_argument(
        "--figure", type=Path, default=Path("figures/rayleigh_taylor_fields.png")
    )
    parser.add_argument(
        "--growth-figure", type=Path, default=Path("figures/rayleigh_taylor_growth.png")
    )
    args = parser.parse_args()
    resolutions = sorted(set(args.resolutions))
    if len(resolutions) < 2 or resolutions[0] < 32:
        raise ValueError("provide at least two resolutions of 32 cells or more")

    histories: list[dict[str, float | int | str]] = []
    summaries: list[dict[str, float | int | str]] = []
    perturbed_states: list[tuple[int, np.ndarray, Grid2D]] = []
    for resolution in resolutions:
        for case in ("control", "perturbed"):
            records, primitive, grid = run_case(
                resolution, case, args.time, args.sample_interval, args.cfl
            )
            histories.extend(records)
            summary = dict(records[-1])
            summary["measured_interface_growth_rate"] = (
                measured_interface_growth_rate(records)
                if case == "perturbed" and args.time >= 2.2
                else ""
            )
            summaries.append(summary)
            if case == "perturbed":
                perturbed_states.append((resolution, primitive, grid))

    args.history.parent.mkdir(parents=True, exist_ok=True)
    with args.history.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(histories[0]))
        writer.writeheader()
        writer.writerows(histories)
    with args.summary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)

    figure, axes = plt.subplots(
        len(perturbed_states), 3, figsize=(12.0, 3.6 * len(perturbed_states)), squeeze=False
    )
    density_minimum = min(float(np.min(state[1][0])) for state in perturbed_states)
    density_maximum = max(float(np.max(state[1][0])) for state in perturbed_states)
    velocity_limit = max(float(np.max(np.abs(state[1][2]))) for state in perturbed_states)
    vorticity_limit = max(
        float(
            np.max(
                np.abs(
                    vorticity_2d(
                        state[1],
                        state[2].dx,
                        state[2].dy,
                        periodic=(True, False),
                    )
                )
            )
        )
        for state in perturbed_states
    )
    for row, (resolution, primitive, grid) in enumerate(perturbed_states):
        vorticity = vorticity_2d(
            primitive, grid.dx, grid.dy, periodic=(True, False)
        )
        fields = (primitive[0], primitive[2], vorticity)
        labels = (r"Density $\rho$", r"Vertical velocity $v_y$", r"Vorticity $\omega_z$")
        color_maps = ("viridis", "RdBu_r", "RdBu_r")
        for column, (field, label, color_map) in enumerate(zip(fields, labels, color_maps)):
            limit = (velocity_limit, vorticity_limit)[column - 1] if column > 0 else None
            image = axes[row, column].pcolormesh(
                grid.x_centers,
                grid.y_centers,
                field.T,
                shading="auto",
                cmap=color_map,
                vmin=-limit if limit is not None else density_minimum,
                vmax=limit if limit is not None else density_maximum,
            )
            axes[row, column].set_aspect("equal")
            axes[row, column].set_xlabel("x")
            axes[row, column].set_ylabel("y")
            axes[row, column].set_title(f"{label}, {resolution}²")
            figure.colorbar(image, ax=axes[row, column], label=label)
    figure.suptitle(f"Rayleigh--Taylor instability at t={args.time:g} (dimensionless)")
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=180, bbox_inches="tight")

    growth_figure, axes = plt.subplots(1, 3, figsize=(14.0, 4.2))
    for resolution in resolutions:
        perturbed = [
            row for row in histories
            if row["resolution"] == resolution and row["case"] == "perturbed"
        ]
        control = [
            row for row in histories
            if row["resolution"] == resolution and row["case"] == "control"
        ]
        times = [float(row["time"]) for row in perturbed]
        axes[0].semilogy(
            times,
            [float(row["transverse_kinetic_energy"]) for row in perturbed],
            label=f"{resolution}²",
        )
        axes[1].plot(times, [float(row["bubble_height"]) for row in perturbed], label=f"bubble {resolution}²")
        axes[1].plot(times, [float(row["spike_depth"]) for row in perturbed], "--", label=f"spike {resolution}²")
        axes[2].semilogy(
            times,
            [float(row["rms_vertical_velocity"]) for row in control],
            label=f"control {resolution}²",
        )
    axes[0].set_ylabel(r"Vertical kinetic energy $K_y$")
    axes[1].set_ylabel("Interface displacement")
    axes[2].set_ylabel(r"Control RMS $v_y$")
    for axis in axes:
        axis.set_xlabel("Time")
        axis.grid(alpha=0.2)
        axis.legend(frameon=False, fontsize=8)
    growth_figure.suptitle("Rayleigh--Taylor growth and hydrostatic drift")
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
