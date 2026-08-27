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


def apply_boundaries_2d(
    state: NDArray[np.float64],
    nghost: int,
    x_boundary: str,
    y_boundary: str,
) -> None:
    """Apply periodic, outflow, or reflective boundaries to a 2D state."""
    if state.ndim != 3 or state.shape[0] != 4:
        raise ValueError("2D state must have shape (4, nx_total, ny_total)")
    if nghost < 1 or min(state.shape[1:]) <= 2 * nghost:
        raise ValueError("2D state must contain active cells and ghost cells")
    valid = ("outflow", "periodic", "reflective")
    if x_boundary not in valid or y_boundary not in valid:
        raise ValueError(f"boundary types must be chosen from {valid}")

    if x_boundary == "outflow":
        state[:, :nghost, :] = state[:, nghost : nghost + 1, :]
        state[:, -nghost:, :] = state[:, -nghost - 1 : -nghost, :]
    elif x_boundary == "periodic":
        state[:, :nghost, :] = state[:, -2 * nghost : -nghost, :]
        state[:, -nghost:, :] = state[:, nghost : 2 * nghost, :]
    else:
        state[:, :nghost, :] = state[:, nghost : 2 * nghost, :][:, ::-1, :]
        state[1, :nghost, :] *= -1.0
        state[:, -nghost:, :] = state[:, -2 * nghost : -nghost, :][:, ::-1, :]
        state[1, -nghost:, :] *= -1.0

    if y_boundary == "outflow":
        state[:, :, :nghost] = state[:, :, nghost : nghost + 1]
        state[:, :, -nghost:] = state[:, :, -nghost - 1 : -nghost]
    elif y_boundary == "periodic":
        state[:, :, :nghost] = state[:, :, -2 * nghost : -nghost]
        state[:, :, -nghost:] = state[:, :, nghost : 2 * nghost]
    else:
        state[:, :, :nghost] = state[:, :, nghost : 2 * nghost][:, :, ::-1]
        state[2, :, :nghost] *= -1.0
        state[:, :, -nghost:] = state[:, :, -2 * nghost : -nghost][:, :, ::-1]
        state[2, :, -nghost:] *= -1.0
