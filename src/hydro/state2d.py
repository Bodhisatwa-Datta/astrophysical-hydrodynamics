"""State conversion and physical fluxes for the two-dimensional Euler system."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .eos import UnphysicalStateError, validate_gamma


def validate_primitive_2d(primitive: ArrayLike) -> None:
    """Reject non-finite or non-positive ``[rho, vx, vy, p]`` states."""
    state = np.asarray(primitive, dtype=np.float64)
    if state.shape[0] != 4:
        raise ValueError("a 2D primitive state must have first dimension of length 4")
    if not np.all(np.isfinite(state)):
        raise UnphysicalStateError("2D primitive state contains NaN or infinity")
    if np.any(state[0] <= 0.0):
        raise UnphysicalStateError("density must be strictly positive")
    if np.any(state[3] <= 0.0):
        raise UnphysicalStateError("pressure must be strictly positive")


def pressure_2d(conserved: ArrayLike, gamma: float) -> NDArray[np.float64]:
    """Return pressure from ``[rho, rho*vx, rho*vy, E]`` without clipping."""
    validate_gamma(gamma)
    state = np.asarray(conserved, dtype=np.float64)
    if state.shape[0] != 4:
        raise ValueError("a 2D conserved state must have first dimension of length 4")
    if not np.all(np.isfinite(state)):
        raise UnphysicalStateError("2D conserved state contains NaN or infinity")
    rho = state[0]
    if np.any(rho <= 0.0):
        raise UnphysicalStateError("density must be strictly positive")
    kinetic = 0.5 * (state[1] ** 2 + state[2] ** 2) / rho
    pressure = (gamma - 1.0) * (state[3] - kinetic)
    if not np.all(np.isfinite(pressure)) or np.any(pressure <= 0.0):
        raise UnphysicalStateError("pressure must be finite and strictly positive")
    return pressure


def primitive_to_conservative_2d(
    primitive: ArrayLike, gamma: float
) -> NDArray[np.float64]:
    """Convert ``[rho, vx, vy, p]`` to ``[rho, rho*vx, rho*vy, E]``."""
    validate_gamma(gamma)
    state = np.asarray(primitive, dtype=np.float64)
    validate_primitive_2d(state)
    rho, velocity_x, velocity_y, pressure = state
    energy = pressure / (gamma - 1.0) + 0.5 * rho * (
        velocity_x**2 + velocity_y**2
    )
    return np.stack(
        (rho, rho * velocity_x, rho * velocity_y, energy), axis=0
    )


def conservative_to_primitive_2d(
    conserved: ArrayLike, gamma: float
) -> NDArray[np.float64]:
    """Convert ``[rho, rho*vx, rho*vy, E]`` to primitive variables."""
    state = np.asarray(conserved, dtype=np.float64)
    pressure = pressure_2d(state, gamma)
    rho = state[0]
    primitive = np.stack(
        (rho, state[1] / rho, state[2] / rho, pressure), axis=0
    )
    validate_primitive_2d(primitive)
    return primitive


def euler_flux_2d(
    conserved: ArrayLike, gamma: float, direction: str
) -> NDArray[np.float64]:
    """Return the physical Euler flux in the x or y direction."""
    state = np.asarray(conserved, dtype=np.float64)
    primitive = conservative_to_primitive_2d(state, gamma)
    rho, velocity_x, velocity_y, pressure = primitive
    if direction == "x":
        return np.stack(
            (
                state[1],
                rho * velocity_x**2 + pressure,
                rho * velocity_x * velocity_y,
                velocity_x * (state[3] + pressure),
            ),
            axis=0,
        )
    if direction == "y":
        return np.stack(
            (
                state[2],
                rho * velocity_x * velocity_y,
                rho * velocity_y**2 + pressure,
                velocity_y * (state[3] + pressure),
            ),
            axis=0,
        )
    raise ValueError("direction must be 'x' or 'y'")

