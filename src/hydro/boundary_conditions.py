"""Ghost-cell boundary conditions."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def apply_outflow(state: NDArray[np.float64], nghost: int) -> None:
    """Fill ghost cells by zero-gradient extrapolation, in place."""
    if state.ndim != 2 or state.shape[0] != 3:
        raise ValueError("state must have shape (3, n_total)")
    if nghost < 1 or state.shape[1] <= 2 * nghost:
        raise ValueError("state must contain at least one active cell and ghost cells")
    state[:, :nghost] = state[:, nghost : nghost + 1]
    state[:, -nghost:] = state[:, -nghost - 1 : -nghost]


def apply_periodic(state: NDArray[np.float64], nghost: int) -> None:
    """Fill ghost cells from the opposite side of the active domain."""
    if state.ndim != 2 or state.shape[0] != 3:
        raise ValueError("state must have shape (3, n_total)")
    if nghost < 1 or state.shape[1] <= 2 * nghost:
        raise ValueError("state must contain at least one active cell and ghost cells")
    state[:, :nghost] = state[:, -2 * nghost : -nghost]
    state[:, -nghost:] = state[:, nghost : 2 * nghost]

