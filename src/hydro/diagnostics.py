"""Reusable diagnostics for one-dimensional hydrodynamic states."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from .eos import sound_speed
from .state import conservative_to_primitive
from .state2d import conservative_to_primitive_2d


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
