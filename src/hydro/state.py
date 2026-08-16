"""Conversions and physical fluxes for the one-dimensional Euler equations."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .eos import pressure, validate_gamma, validate_primitive


def primitive_to_conservative(primitive: ArrayLike, gamma: float) -> NDArray[np.float64]:
    """Convert ``[rho, u, p]`` to ``[rho, rho*u, E]``."""
    validate_gamma(gamma)
    state = np.asarray(primitive, dtype=np.float64)
    if state.shape[0] != 3:
        raise ValueError("a 1D primitive state must have first dimension of length 3")
    rho, velocity, p = state
    validate_primitive(rho, velocity, p)
    energy = p / (gamma - 1.0) + 0.5 * rho * velocity**2
    return np.stack((rho, rho * velocity, energy), axis=0)


def conservative_to_primitive(conserved: ArrayLike, gamma: float) -> NDArray[np.float64]:
    """Convert ``[rho, rho*u, E]`` to ``[rho, u, p]``."""
    state = np.asarray(conserved, dtype=np.float64)
    if state.shape[0] != 3:
        raise ValueError("a 1D conserved state must have first dimension of length 3")
    p = pressure(state, gamma)
    rho = state[0]
    velocity = state[1] / rho
    validate_primitive(rho, velocity, p)
    return np.stack((rho, velocity, p), axis=0)


def euler_flux(conserved: ArrayLike, gamma: float) -> NDArray[np.float64]:
    """Return ``[rho*u, rho*u^2+p, u*(E+p)]`` for a valid state."""
    state = np.asarray(conserved, dtype=np.float64)
    primitive = conservative_to_primitive(state, gamma)
    rho, velocity, p = primitive
    return np.stack(
        (state[1], rho * velocity**2 + p, velocity * (state[2] + p)), axis=0
    )
