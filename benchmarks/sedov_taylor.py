"""Run a 2D Sedov--Taylor blast resolution and similarity study."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from time import perf_counter

import matplotlib
import numpy as np

from hydro.diagnostics import radial_profile_2d, shock_radius_2d, state_summary_2d, totals_2d
from hydro.problems import sedov_taylor_2d
from hydro.sedov import sedov_similarity_2d
from hydro.solver2d import Grid2D, Solver2D

matplotlib.use("Agg")
import matplotlib.pyplot as plt


GAMMA = 1.4
EXPLOSION_ENERGY = 1.0
INJECTION_RADIUS = 0.05


def sector_shock_radii(
    primitive: np.ndarray,
    grid: Grid2D,
    sectors: int = 16,
) -> np.ndarray:
    """Measure the density-gradient shock radius in equal angular sectors."""
    x, y = grid.mesh
    radius = np.sqrt(x**2 + y**2)
    angle = (np.arctan2(y, x) + 2.0 * np.pi) % (2.0 * np.pi)
    radii = np.empty(sectors)
    bin_width = min(grid.dx, grid.dy)
    for sector in range(sectors):
        lower = 2.0 * np.pi * sector / sectors
        upper = 2.0 * np.pi * (sector + 1) / sectors
        mask = (angle >= lower) & (angle < upper) & (radius <= 0.45)
        radial_coordinate = radius[mask]
        density = primitive[0][mask]
        n_bins = int(np.ceil(0.45 / bin_width))
        indices = np.floor(radial_coordinate / bin_width).astype(int)
        counts = np.bincount(indices, minlength=n_bins)
        sums = np.bincount(indices, weights=density, minlength=n_bins)
        populated = counts > 0
        bin_centers = (np.arange(n_bins) + 0.5) * bin_width
        profile = np.full(n_bins, np.nan)
        profile[populated] = sums[populated] / counts[populated]
        valid = populated & (bin_centers >= INJECTION_RADIUS)
        valid_centers = bin_centers[valid]
        gradient = np.gradient(profile[valid], valid_centers)
        local_index = int(np.argmin(gradient))
        shock_radius = float(valid_centers[local_index])
        if 0 < local_index < gradient.size - 1:
            left, middle, right = gradient[local_index - 1 : local_index + 2]
            denominator = left - 2.0 * middle + right
            if denominator != 0.0:
                offset = 0.5 * (left - right) / denominator
                if abs(offset) <= 1.0:
                    shock_radius += float(offset * bin_width)
        radii[sector] = shock_radius
    return radii


def record(
    solver: Solver2D,
    initial_totals: dict[str, float],
    runtime: float,
) -> dict[str, float | int]:
    grid = solver.grid
    totals = totals_2d(solver.active_conserved, grid.dx * grid.dy)
    primitive = solver.primitive
    measured_shock_radius = shock_radius_2d(
        primitive,
        *grid.mesh,
        bin_width=min(grid.dx, grid.dy),
        minimum_radius=INJECTION_RADIUS,
        maximum_radius=0.45,
    )
    exact_shock_radius = sedov_similarity_2d(
        np.asarray([0.0]),
        solver.time,
        gamma=solver.gamma,
        explosion_energy=EXPLOSION_ENERGY,
        ambient_pressure=1.0e-5,
    ).shock_radius
    diagnostics: dict[str, float | int] = {
        "resolution": grid.nx,
        "time": solver.time,
        "steps": solver.steps,
        "runtime_seconds": runtime,
        "shock_radius": measured_shock_radius,
        "exact_shock_radius": exact_shock_radius,
        "relative_shock_radius_error": (
            measured_shock_radius - exact_shock_radius
        )
        / exact_shock_radius,
        "relative_mass_change": (
            totals["mass"] - initial_totals["mass"]
        )
        / initial_totals["mass"],
        "relative_energy_change": (
            totals["energy"] - initial_totals["energy"]
        )
        / initial_totals["energy"],
        "absolute_momentum_x_change": totals["momentum_x"]
        - initial_totals["momentum_x"],
        "absolute_momentum_y_change": totals["momentum_y"]
        - initial_totals["momentum_y"],
    }
    diagnostics.update(state_summary_2d(solver.active_conserved, solver.gamma))
    return diagnostics


def run_resolution(
    resolution: int,
    final_time: float,
    sample_interval: float,
    cfl: float,
) -> tuple[list[dict[str, float | int]], np.ndarray, Grid2D]:
    grid = Grid2D(-0.5, 0.5, resolution, -0.5, 0.5, resolution)
    solver = Solver2D(
        grid,
        gamma=GAMMA,
        cfl=cfl,
        riemann_solver="hll",
        reconstruction="muscl",
        limiter="mc",
        integrator="rk2",
        x_boundary="outflow",
        y_boundary="outflow",
    )
    solver.initialise(
        sedov_taylor_2d(
            *grid.mesh,
            cell_area=grid.dx * grid.dy,
            gamma=GAMMA,
            explosion_energy=EXPLOSION_ENERGY,
            injection_radius=INJECTION_RADIUS,
        )
    )
    initial_totals = totals_2d(solver.active_conserved, grid.dx * grid.dy)
    records: list[dict[str, float | int]] = []
    start = perf_counter()
    sample_index = 1
    target = min(sample_index * sample_interval, final_time)
    while solver.time < final_time:
        solver.step(min(solver.timestep(), target - solver.time))
        if np.isclose(solver.time, target, rtol=0.0, atol=32.0 * np.finfo(float).eps):
            records.append(record(solver, initial_totals, perf_counter() - start))
            if np.isclose(target, final_time, rtol=0.0, atol=32.0 * np.finfo(float).eps):
                break
            sample_index += 1
            target = min(sample_index * sample_interval, final_time)
    return records, solver.primitive.copy(), grid


def fit_similarity(records: list[dict[str, float | int]]) -> tuple[float, float]:
    """Fit ``R=C*t**alpha`` after the finite injection region is forgotten."""
    selected = [row for row in records if float(row["time"]) >= 0.015]
    times = np.asarray([float(row["time"]) for row in selected])
    radii = np.asarray([float(row["shock_radius"]) for row in selected])
    exponent, log_coefficient = np.polyfit(np.log(times), np.log(radii), 1)
    return float(exponent), float(np.exp(log_coefficient))


def similarity_profile_errors(
    primitive: np.ndarray,
    grid: Grid2D,
    time: float,
) -> dict[str, float]:
    """Compare radial means with the exact strong-shock similarity profiles."""
    x, y = grid.mesh
    radius_2d = np.sqrt(x**2 + y**2)
    radial_velocity_2d = np.divide(
        x * primitive[1] + y * primitive[2],
        radius_2d,
        out=np.zeros_like(radius_2d),
        where=radius_2d > 0.0,
    )
    profiles = []
    radii = np.empty(0)
    for field in (primitive[0], radial_velocity_2d, primitive[3]):
        radii, profile, counts = radial_profile_2d(
            field,
            x,
            y,
            bin_width=min(grid.dx, grid.dy),
            maximum_radius=0.45,
        )
        profiles.append(profile)
    exact = sedov_similarity_2d(
        radii,
        time,
        gamma=GAMMA,
        explosion_energy=EXPLOSION_ENERGY,
        ambient_pressure=1.0e-5,
    )
    comparison = (radii <= 1.25 * exact.shock_radius) & (counts > 0)
    exact_profiles = (exact.density, exact.radial_velocity, exact.pressure)
    names = ("density", "radial_velocity", "pressure")
    errors = {}
    for name, numerical, reference in zip(names, profiles, exact_profiles):
        compared_radii = radii[comparison]
        absolute_error = (
            np.abs(numerical[comparison] - reference[comparison]) * compared_radii
        )
        reference_magnitude = np.abs(reference[comparison]) * compared_radii
        radial_widths = np.diff(compared_radii)
        numerator = np.sum(
            0.5 * (absolute_error[:-1] + absolute_error[1:]) * radial_widths
        )
        denominator = np.sum(
            0.5
            * (reference_magnitude[:-1] + reference_magnitude[1:])
            * radial_widths
        )
        errors[f"{name}_profile_relative_l1"] = float(numerator / denominator)
    return errors


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[64, 96, 128])
    parser.add_argument("--time", type=float, default=0.05)
    parser.add_argument("--sample-interval", type=float, default=0.005)
    parser.add_argument("--cfl", type=float, default=0.25)
    parser.add_argument(
        "--history", type=Path, default=Path("benchmarks/sedov_taylor_history.csv")
    )
    parser.add_argument(
        "--summary", type=Path, default=Path("benchmarks/sedov_taylor_summary.csv")
    )
    parser.add_argument(
        "--figure", type=Path, default=Path("figures/sedov_taylor_fields.png")
    )
    parser.add_argument(
        "--analysis-figure", type=Path, default=Path("figures/sedov_taylor_analysis.png")
    )
    args = parser.parse_args()
    resolutions = sorted(set(args.resolutions))
    if len(resolutions) < 2 or resolutions[0] < 32:
        raise ValueError("provide at least two resolutions of 32 cells or more")

    histories: list[dict[str, float | int]] = []
    summaries: list[dict[str, float | int]] = []
    final_states: list[tuple[int, np.ndarray, Grid2D]] = []
    for resolution in resolutions:
        records, primitive, grid = run_resolution(
            resolution, args.time, args.sample_interval, args.cfl
        )
        histories.extend(records)
        exponent, coefficient = fit_similarity(records)
        exact_coefficient = sedov_similarity_2d(
            np.asarray([0.0]),
            args.time,
            gamma=GAMMA,
            explosion_energy=EXPLOSION_ENERGY,
            ambient_pressure=1.0e-5,
        ).shock_radius / np.sqrt(args.time)
        angular_radii = sector_shock_radii(primitive, grid)
        summary = dict(records[-1])
        summary.update(
            {
                "similarity_exponent": exponent,
                "similarity_coefficient": coefficient,
                "exact_similarity_coefficient": exact_coefficient,
                "relative_similarity_coefficient_error": (
                    coefficient - exact_coefficient
                )
                / exact_coefficient,
                "relative_exponent_error": (exponent - 0.5) / 0.5,
                "angular_shock_radius_mean": float(np.mean(angular_radii)),
                "angular_shock_radius_relative_std": float(
                    np.std(angular_radii) / np.mean(angular_radii)
                ),
                "angular_shock_radius_peak_to_peak": float(
                    (np.max(angular_radii) - np.min(angular_radii))
                    / np.mean(angular_radii)
                ),
            }
        )
        summary.update(similarity_profile_errors(primitive, grid, args.time))
        summaries.append(summary)
        final_states.append((resolution, primitive, grid))

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
        len(final_states), 2, figsize=(8.5, 4.0 * len(final_states)), squeeze=False
    )
    density_min = min(float(np.min(item[1][0])) for item in final_states)
    density_max = max(float(np.max(item[1][0])) for item in final_states)
    pressure_min = min(float(np.min(item[1][3])) for item in final_states)
    pressure_max = max(float(np.max(item[1][3])) for item in final_states)
    for row, (resolution, primitive, grid) in enumerate(final_states):
        for column, (field, label, limits, color_map) in enumerate(
            (
                (primitive[0], r"Density $\rho$", (density_min, density_max), "viridis"),
                (primitive[3], "Pressure p", (pressure_min, pressure_max), "magma"),
            )
        ):
            image = axes[row, column].pcolormesh(
                grid.x_centers,
                grid.y_centers,
                field.T,
                shading="auto",
                cmap=color_map,
                vmin=limits[0],
                vmax=limits[1],
            )
            axes[row, column].set_aspect("equal")
            axes[row, column].set_xlabel("x")
            axes[row, column].set_ylabel("y")
            axes[row, column].set_title(f"{label}, {resolution}²")
            figure.colorbar(image, ax=axes[row, column], label=label)
    figure.suptitle(f"Two-dimensional Sedov--Taylor blast at t={args.time:g}")
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=180, bbox_inches="tight")

    analysis, axes = plt.subplots(2, 2, figsize=(11.5, 8.5))
    axes = axes.ravel()
    reference_times = np.linspace(0.015, args.time, 100)
    for (resolution, primitive, grid), summary in zip(final_states, summaries):
        rows = [row for row in histories if row["resolution"] == resolution]
        axes[0].plot(
            [float(row["time"]) for row in rows],
            [float(row["shock_radius"]) for row in rows],
            "o-",
            label=f"{resolution}²",
        )
        coefficient = float(summary["similarity_coefficient"])
        exponent = float(summary["similarity_exponent"])
        axes[0].plot(
            reference_times,
            coefficient * reference_times**exponent,
            "--",
            alpha=0.55,
            label="_nolegend_",
        )
        radii, density, _ = radial_profile_2d(
            primitive[0], *grid.mesh, bin_width=min(grid.dx, grid.dy), maximum_radius=0.45
        )
        _, pressure, _ = radial_profile_2d(
            primitive[3], *grid.mesh, bin_width=min(grid.dx, grid.dy), maximum_radius=0.45
        )
        x, y = grid.mesh
        radius_2d = np.sqrt(x**2 + y**2)
        radial_velocity_2d = np.divide(
            x * primitive[1] + y * primitive[2],
            radius_2d,
            out=np.zeros_like(radius_2d),
            where=radius_2d > 0.0,
        )
        _, radial_velocity, _ = radial_profile_2d(
            radial_velocity_2d,
            x,
            y,
            bin_width=min(grid.dx, grid.dy),
            maximum_radius=0.45,
        )
        axes[1].plot(radii, density, label=f"{resolution}²")
        axes[2].semilogy(radii, pressure, label=f"{resolution}²")
        axes[3].plot(radii, radial_velocity, label=f"{resolution}²")
    reference = sedov_similarity_2d(
        np.linspace(0.0, 0.45, 1200),
        args.time,
        gamma=GAMMA,
        explosion_energy=EXPLOSION_ENERGY,
        ambient_pressure=1.0e-5,
    )
    axes[0].plot(
        reference_times,
        reference.shock_radius * np.sqrt(reference_times / args.time),
        "k:",
        linewidth=2.0,
        label="exact similarity solution",
    )
    axes[1].plot(reference.radius, reference.density, "k:", linewidth=2.0, label="exact")
    axes[2].semilogy(
        reference.radius, reference.pressure, "k:", linewidth=2.0, label="exact"
    )
    axes[3].plot(
        reference.radius,
        reference.radial_velocity,
        "k:",
        linewidth=2.0,
        label="exact",
    )
    axes[0].set_ylabel("Shock radius")
    axes[0].set_xlabel("Time")
    axes[0].set_title("Measured radius and power-law fits")
    axes[1].set_xlabel("Radius")
    axes[1].set_ylabel(r"Mean density $\rho$")
    axes[1].set_title(f"Radial density at t={args.time:g}")
    axes[2].set_xlabel("Radius")
    axes[2].set_ylabel("Mean pressure p")
    axes[2].set_title(f"Radial pressure at t={args.time:g}")
    axes[3].set_xlabel("Radius")
    axes[3].set_ylabel("Mean radial velocity")
    axes[3].set_title(f"Radial velocity at t={args.time:g}")
    for axis in axes:
        axis.grid(alpha=0.2)
        axis.legend(frameon=False)
    analysis.tight_layout()
    args.analysis_figure.parent.mkdir(parents=True, exist_ok=True)
    analysis.savefig(args.analysis_figure, dpi=180, bbox_inches="tight")

    print(json.dumps(summaries, indent=2))
    print(f"Saved {args.history}")
    print(f"Saved {args.summary}")
    print(f"Saved {args.figure}")
    print(f"Saved {args.analysis_figure}")


if __name__ == "__main__":
    main()
