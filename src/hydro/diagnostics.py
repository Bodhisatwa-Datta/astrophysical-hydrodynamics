"""Reusable diagnostics for one-dimensional hydrodynamic states."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from .eos import sound_speed
from .state import conservative_to_primitive


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

