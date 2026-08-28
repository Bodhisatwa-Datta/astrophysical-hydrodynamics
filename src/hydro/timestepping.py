"""CFL stability limits for explicit finite-volume integration."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from .eos import sound_speed
from .state import conservative_to_primitive


def cfl_timestep(
    conserved: ArrayLike, dx: float, gamma: float, cfl: float, nghost: int = 0
) -> float:
    """Return ``C_CFL*dx/max(|u|+c_s)`` over active cells."""
    if dx <= 0.0 or not np.isfinite(dx):
        raise ValueError("dx must be finite and positive")
    if not 0.0 < cfl <= 1.0:
        raise ValueError("cfl must lie in (0, 1]")
    state = np.asarray(conserved, dtype=np.float64)
    active = state[:, nghost:-nghost] if nghost else state
    primitive = conservative_to_primitive(active, gamma)
    speed = np.abs(primitive[1]) + sound_speed(primitive[0], primitive[2], gamma)
    maximum_speed = float(np.max(speed))
    if maximum_speed <= 0.0 or not np.isfinite(maximum_speed):
        raise ValueError("maximum signal speed must be finite and positive")
    return cfl * dx / maximum_speed
