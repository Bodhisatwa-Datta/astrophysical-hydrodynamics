"""Ideal-gas equation of state and physical-state validation."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


class UnphysicalStateError(ValueError):
    """Raised when a hydrodynamic state is non-finite or non-positive."""


def _as_float_array(value: ArrayLike) -> NDArray[np.float64]:
    return np.asarray(value, dtype=np.float64)


def validate_gamma(gamma: float) -> None:
    """Validate the ideal-gas adiabatic index."""
    if gamma <= 1.0 or not np.isfinite(gamma):
        raise ValueError("gamma must be finite and greater than one")


def validate_primitive(rho: ArrayLike, velocity: ArrayLike, p: ArrayLike) -> None:
    """Raise :class:`UnphysicalStateError` for an invalid primitive state."""
    rho_arr, velocity_arr, p_arr = map(_as_float_array, (rho, velocity, p))
    if not (
        np.all(np.isfinite(rho_arr))
        and np.all(np.isfinite(velocity_arr))
        and np.all(np.isfinite(p_arr))
    ):
        raise UnphysicalStateError("primitive state contains NaN or infinity")
    if np.any(rho_arr <= 0.0):
        raise UnphysicalStateError("density must be strictly positive")
    if np.any(p_arr <= 0.0):
        raise UnphysicalStateError("pressure must be strictly positive")


def pressure(conserved: ArrayLike, gamma: float) -> NDArray[np.float64]:
    """Return pressure from ``[rho, rho*u, E]`` without clipping."""
    validate_gamma(gamma)
    state = _as_float_array(conserved)
    if state.shape[0] != 3:
        raise ValueError("a 1D conserved state must have first dimension of length 3")
    if not np.all(np.isfinite(state)):
        raise UnphysicalStateError("conserved state contains NaN or infinity")
    rho = state[0]
    if np.any(rho <= 0.0):
        raise UnphysicalStateError("density must be strictly positive")
    kinetic = 0.5 * state[1] ** 2 / rho
    p = (gamma - 1.0) * (state[2] - kinetic)
    if not np.all(np.isfinite(p)) or np.any(p <= 0.0):
        raise UnphysicalStateError("pressure must be finite and strictly positive")
    return p


def sound_speed(rho: ArrayLike, p: ArrayLike, gamma: float) -> NDArray[np.float64]:
    """Return the ideal-gas adiabatic sound speed ``sqrt(gamma*p/rho)``."""
    validate_gamma(gamma)
    rho_arr, p_arr = map(_as_float_array, (rho, p))
    validate_primitive(rho_arr, np.zeros_like(rho_arr + p_arr), p_arr)
    return np.sqrt(gamma * p_arr / rho_arr)
