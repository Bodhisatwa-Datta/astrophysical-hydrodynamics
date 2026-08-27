"""Reusable execution and error analysis for 1D Riemann benchmarks."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

import numpy as np
from numpy.typing import NDArray

from .diagnostics import state_summary, totals
from .exact_riemann import exact_riemann_solution
from .problems import RiemannProblem, entropy_wave
from .solver1d import Grid1D, Solver1D
from .state import euler_flux, primitive_to_conservative


@dataclass(frozen=True)
class RiemannResult:
    """Numerical and exact states plus measured diagnostics."""

    problem: RiemannProblem
    x: NDArray[np.float64]
    numerical: NDArray[np.float64]
    exact: NDArray[np.float64]
    diagnostics: dict[str, float | int | str]


@dataclass(frozen=True)
class SmoothWaveResult:
    """Numerical and exact states for one periodic entropy-wave crossing."""

    x: NDArray[np.float64]
    numerical: NDArray[np.float64]
    exact: NDArray[np.float64]
    diagnostics: dict[str, float | int | str]


def run_riemann_problem(
    problem: RiemannProblem,
    n_cells: int,
    cfl: float = 0.8,
    reconstruction: str = "constant",
    limiter: str = "mc",
    integrator: str = "euler",
) -> RiemannResult:
    """Run one problem on ``[0, 1]`` and compare with its exact solution."""
    grid = Grid1D(0.0, 1.0, n_cells)
    solver = Solver1D(
        grid,
        gamma=problem.gamma,
        cfl=cfl,
        reconstruction=reconstruction,
        limiter=limiter,
        integrator=integrator,
    )
    solver.initialise_function(problem.initial_condition)
    initial_totals = totals(solver.active_conserved, grid.dx)
    start = perf_counter()
    solver.run(problem.final_time)
    runtime = perf_counter() - start
    final_totals = totals(solver.active_conserved, grid.dx)
    left_flux = euler_flux(
        primitive_to_conservative(np.asarray(problem.left), problem.gamma), problem.gamma
    )
    right_flux = euler_flux(
        primitive_to_conservative(np.asarray(problem.right), problem.gamma), problem.gamma
    )
    expected_boundary_change = (left_flux - right_flux) * problem.final_time
    numerical = solver.primitive.copy()
    exact = exact_riemann_solution(
        grid.centers,
        problem.final_time,
        problem.left,
        problem.right,
        problem.gamma,
        problem.discontinuity,
    )
    field_names = ("density", "velocity", "pressure")
    diagnostics: dict[str, float | int | str] = {
        "cells": n_cells,
        "time": solver.time,
        "steps": solver.steps,
        "runtime_seconds": runtime,
        "reconstruction": reconstruction,
        "limiter": limiter if reconstruction == "muscl" else "none",
        "integrator": integrator,
    }
    for index, name in enumerate(field_names):
        difference = numerical[index] - exact[index]
        diagnostics[f"{name}_L1_error"] = float(grid.dx * np.sum(np.abs(difference)))
        diagnostics[f"{name}_Linf_error"] = float(np.max(np.abs(difference)))
    for index, quantity in enumerate(("mass", "momentum", "energy")):
        initial = initial_totals[quantity]
        change = final_totals[quantity] - initial
        residual = change - float(expected_boundary_change[index])
        diagnostics[f"absolute_{quantity}_change"] = change
        diagnostics[f"{quantity}_boundary_budget_residual"] = residual
        if initial != 0.0:
            diagnostics[f"relative_{quantity}_change"] = change / abs(initial)
            diagnostics[f"relative_{quantity}_boundary_budget_residual"] = (
                residual / abs(initial)
            )
    diagnostics.update(state_summary(solver.active_conserved, problem.gamma))
    return RiemannResult(problem, grid.centers, numerical, exact, diagnostics)


def contact_transition_width(result: RiemannResult, fraction: float = 0.01) -> float:
    """Width containing the central 98% of a density contact transition."""
    if not 0.0 < fraction < 0.5:
        raise ValueError("fraction must lie in (0, 0.5)")
    rho_left = result.problem.left[0]
    rho_right = result.problem.right[0]
    low, high = sorted((rho_left, rho_right))
    span = high - low
    transition = (result.numerical[0] > low + fraction * span) & (
        result.numerical[0] < high - fraction * span
    )
    indices = np.flatnonzero(transition)
    if indices.size == 0:
        return 0.0
    dx = float(result.x[1] - result.x[0])
    return float((indices[-1] - indices[0] + 1) * dx)


def run_entropy_wave(
    n_cells: int,
    cfl: float,
    reconstruction: str,
    limiter: str,
    integrator: str,
    final_time: float = 1.0,
) -> SmoothWaveResult:
    """Advect a smooth entropy wave on a periodic unit domain."""
    grid = Grid1D(0.0, 1.0, n_cells)
    solver = Solver1D(
        grid,
        gamma=1.4,
        cfl=cfl,
        reconstruction=reconstruction,
        limiter=limiter,
        integrator=integrator,
        boundary="periodic",
    )
    solver.initialise_function(entropy_wave)
    initial_totals = totals(solver.active_conserved, grid.dx)
    start = perf_counter()
    solver.run(final_time)
    runtime = perf_counter() - start
    final_totals = totals(solver.active_conserved, grid.dx)
    numerical = solver.primitive.copy()
    exact = entropy_wave(grid.centers, time=final_time)
    difference = numerical - exact
    diagnostics: dict[str, float | int | str] = {
        "cells": n_cells,
        "time": solver.time,
        "steps": solver.steps,
        "runtime_seconds": runtime,
        "reconstruction": reconstruction,
        "limiter": limiter if reconstruction == "muscl" else "none",
        "integrator": integrator,
        "density_L1_error": float(grid.dx * np.sum(np.abs(difference[0]))),
        "density_Linf_error": float(np.max(np.abs(difference[0]))),
    }
    for quantity in ("mass", "momentum", "energy"):
        initial = initial_totals[quantity]
        diagnostics[f"relative_{quantity}_change"] = (
            final_totals[quantity] - initial
        ) / abs(initial)
    diagnostics.update(state_summary(solver.active_conserved, 1.4))
    return SmoothWaveResult(grid.centers, numerical, exact, diagnostics)
