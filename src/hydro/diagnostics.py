"""Reusable diagnostics for one-dimensional hydrodynamic states."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from .eos import sound_speed
from .state import conservative_to_primitive
from .state2d import conservative_to_primitive_2d
from .gravity import gravitational_potential_2d


def totals(conserved: ArrayLike, cell_volume: float) -> dict[str, float]:
    """Return finite-volume totals of mass, momentum, and total energy."""
    state = np.asarray(conserved, dtype=np.float64)
    if state.shape[0] != 3 or cell_volume <= 0.0:
        raise ValueError("expected a (3, ...) state and positive cell volume")
    values = np.sum(state, axis=tuple(range(1, state.ndim))) * cell_volume
    return {"mass": float(values[0]), "momentum": float(values[1]), "energy": float(values[2])}


def state_summary(conserved: ArrayLike, gamma: float) -> dict[str, float]:
    """Return minimum thermodynamic values and maximum Mach number."""
    primitive = conservative_to_primitive(conserved, gamma)
    rho, velocity, p = primitive
    mach = np.abs(velocity) / sound_speed(rho, p, gamma)
    return {
        "minimum_density": float(np.min(rho)),
        "minimum_pressure": float(np.min(p)),
        "maximum_mach": float(np.max(mach)),
    }


def totals_2d(conserved: ArrayLike, cell_area: float) -> dict[str, float]:
    """Return 2D totals of mass, both momenta, and total energy."""
    state = np.asarray(conserved, dtype=np.float64)
    if state.shape[0] != 4 or cell_area <= 0.0:
        raise ValueError("expected a (4, ...) state and positive cell area")
    values = np.sum(state, axis=tuple(range(1, state.ndim))) * cell_area
    return {
        "mass": float(values[0]),
        "momentum_x": float(values[1]),
        "momentum_y": float(values[2]),
        "energy": float(values[3]),
    }


def state_summary_2d(conserved: ArrayLike, gamma: float) -> dict[str, float]:
    """Return minimum thermodynamic values and maximum 2D Mach number."""
    primitive = conservative_to_primitive_2d(conserved, gamma)
    density, velocity_x, velocity_y, pressure = primitive
    speed = np.sqrt(velocity_x**2 + velocity_y**2)
    sound = sound_speed(density, pressure, gamma)
    return {
        "minimum_density": float(np.min(density)),
        "minimum_pressure": float(np.min(pressure)),
        "maximum_mach": float(np.max(speed / sound)),
    }


def vorticity_2d(
    primitive: ArrayLike,
    dx: float,
    dy: float,
    periodic: bool | tuple[bool, bool] = True,
) -> np.ndarray:
    """Return cell-centred z-vorticity, ``d(vy)/dx - d(vx)/dy``."""
    state = np.asarray(primitive, dtype=np.float64)
    if state.ndim != 3 or state.shape[0] != 4:
        raise ValueError("expected a primitive state with shape (4, nx, ny)")
    if dx <= 0.0 or dy <= 0.0:
        raise ValueError("cell spacings must be positive")
    velocity_x = state[1]
    velocity_y = state[2]
    periodic_x, periodic_y = (periodic, periodic) if isinstance(periodic, bool) else periodic
    if periodic_x:
        dvy_dx = (np.roll(velocity_y, -1, axis=0) - np.roll(velocity_y, 1, axis=0)) / (
            2.0 * dx
        )
    else:
        dvy_dx = np.gradient(velocity_y, dx, axis=0, edge_order=2)
    if periodic_y:
        dvx_dy = (np.roll(velocity_x, -1, axis=1) - np.roll(velocity_x, 1, axis=1)) / (
            2.0 * dy
        )
    else:
        dvx_dy = np.gradient(velocity_x, dy, axis=1, edge_order=2)
    return dvy_dx - dvx_dy


def transverse_kinetic_energy_2d(
    primitive: ArrayLike, cell_area: float
) -> float:
    """Return the integrated kinetic energy in the y-velocity component."""
    state = np.asarray(primitive, dtype=np.float64)
    if state.ndim != 3 or state.shape[0] != 4 or cell_area <= 0.0:
        raise ValueError("expected a (4, nx, ny) state and positive cell area")
    return float(0.5 * cell_area * np.sum(state[0] * state[2] ** 2))


def total_energy_with_gravity_2d(
    conserved: ArrayLike,
    x: ArrayLike,
    y: ArrayLike,
    gravity: tuple[float, float],
    cell_area: float,
) -> float:
    """Return integrated gas energy plus uniform-field potential energy."""
    state = np.asarray(conserved, dtype=np.float64)
    potential = gravitational_potential_2d(x, y, gravity)
    if state.ndim != 3 or state.shape[0] != 4:
        raise ValueError("expected a conserved state with shape (4, nx, ny)")
    if state.shape[1:] != potential.shape or cell_area <= 0.0:
        raise ValueError("state and coordinate shapes must match with positive area")
    return float(cell_area * np.sum(state[3] + state[0] * potential))


def density_interface_height_2d(
    primitive: ArrayLike,
    y_centers: ArrayLike,
    threshold: float = 1.5,
) -> np.ndarray:
    """Interpolate the first upward density-threshold crossing in each column."""
    state = np.asarray(primitive, dtype=np.float64)
    coordinates = np.asarray(y_centers, dtype=np.float64)
    if state.ndim != 3 or state.shape[0] != 4:
        raise ValueError("expected a primitive state with shape (4, nx, ny)")
    if coordinates.ndim != 1 or coordinates.size != state.shape[2]:
        raise ValueError("y_centers must match the state's y dimension")
    heights = np.empty(state.shape[1], dtype=np.float64)
    for index, column in enumerate(state[0]):
        crossings = np.flatnonzero((column[:-1] < threshold) & (column[1:] >= threshold))
        if crossings.size == 0:
            heights[index] = np.nan
            continue
        lower = int(crossings[0])
        fraction = (threshold - column[lower]) / (column[lower + 1] - column[lower])
        heights[index] = coordinates[lower] + fraction * (
            coordinates[lower + 1] - coordinates[lower]
        )
    return heights


def radial_profile_2d(
    field: ArrayLike,
    x: ArrayLike,
    y: ArrayLike,
    bin_width: float,
    center: tuple[float, float] = (0.0, 0.0),
    maximum_radius: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return annular cell averages, bin centres, and cell counts."""
    values = np.asarray(field, dtype=np.float64)
    x_coordinates = np.asarray(x, dtype=np.float64)
    y_coordinates = np.asarray(y, dtype=np.float64)
    if values.shape != x_coordinates.shape or values.shape != y_coordinates.shape:
        raise ValueError("field and coordinate arrays must have matching shapes")
    if bin_width <= 0.0:
        raise ValueError("bin_width must be positive")
    radius = np.sqrt(
        (x_coordinates - center[0]) ** 2 + (y_coordinates - center[1]) ** 2
    )
    if maximum_radius is None:
        maximum_radius = float(np.max(radius)) + 0.5 * bin_width
    n_bins = max(1, int(np.ceil(maximum_radius / bin_width)))
    indices = np.floor(radius / bin_width).astype(int)
    mask = indices < n_bins
    counts = np.bincount(indices[mask].ravel(), minlength=n_bins)
    sums = np.bincount(
        indices[mask].ravel(), weights=values[mask].ravel(), minlength=n_bins
    )
    averages = np.full(n_bins, np.nan, dtype=np.float64)
    populated = counts > 0
    averages[populated] = sums[populated] / counts[populated]
    centers = (np.arange(n_bins, dtype=np.float64) + 0.5) * bin_width
    return centers, averages, counts


def shock_radius_2d(
    primitive: ArrayLike,
    x: ArrayLike,
    y: ArrayLike,
    bin_width: float,
    minimum_radius: float = 0.05,
    maximum_radius: float | None = None,
) -> float:
    """Locate the blast shock from the steepest outward density decrease."""
    state = np.asarray(primitive, dtype=np.float64)
    if state.ndim != 3 or state.shape[0] != 4:
        raise ValueError("expected a primitive state with shape (4, nx, ny)")
    radii, density, counts = radial_profile_2d(
        state[0], x, y, bin_width, maximum_radius=maximum_radius
    )
    valid = np.isfinite(density) & (counts > 0)
    if np.count_nonzero(valid) < 3:
        raise ValueError("too few populated radial bins to locate a shock")
    gradient = np.gradient(density[valid], radii[valid])
    valid_radii = radii[valid]
    candidates = np.flatnonzero(valid_radii >= minimum_radius)
    if candidates.size == 0:
        raise ValueError("no radial bins lie beyond minimum_radius")
    index = int(candidates[np.argmin(gradient[candidates])])
    radius = float(valid_radii[index])
    if 0 < index < gradient.size - 1:
        left, middle, right = gradient[index - 1 : index + 2]
        denominator = left - 2.0 * middle + right
        if denominator != 0.0:
            offset = 0.5 * (left - right) / denominator
            if abs(offset) <= 1.0:
                radius += float(offset * bin_width)
    return radius
