"""Standard initial conditions for one-dimensional validation problems."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def sod_initial_condition(
    x: NDArray[np.float64], discontinuity: float = 0.5
) -> NDArray[np.float64]:
    """Return the standard Sod states ``(1,0,1)`` and ``(0.125,0,0.1)``."""
    left = x < discontinuity
    return np.stack(
        (
            np.where(left, 1.0, 0.125),
            np.zeros_like(x),
            np.where(left, 1.0, 0.1),
        ),
        axis=0,
    )

