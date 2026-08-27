"""Uniform-gravity source terms and potential-energy helpers."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def validate_gravity(gravity: tuple[float, float]) -> tuple[float, float]:
    """Return a finite two-component gravitational acceleration."""
    if len(gravity) != 2:
        raise ValueError("gravity must contain x and y accelerations")
    acceleration = tuple(float(value) for value in gravity)
    if not all(np.isfinite(value) for value in acceleration):
        raise ValueError("gravity components must be finite")
    return acceleration


def constant_gravity_source_2d(
    conserved: ArrayLike, gravity: tuple[float, float]
) -> NDArray[np.float64]:
    """Return ``[0, rho gx, rho gy, momentum dot gravity]``."""
    state = np.asarray(conserved, dtype=np.float64)
    if state.ndim < 1 or state.shape[0] != 4:
        raise ValueError("expected a conserved state with first dimension 4")
    gravity_x, gravity_y = validate_gravity(gravity)
    source = np.zeros_like(state)
    source[1] = state[0] * gravity_x
    source[2] = state[0] * gravity_y
    source[3] = state[1] * gravity_x + state[2] * gravity_y
    return source


def gravitational_potential_2d(
    x: ArrayLike, y: ArrayLike, gravity: tuple[float, float]
) -> NDArray[np.float64]:
    """Return the potential ``Phi=-(gx*x + gy*y)`` for uniform gravity."""
    x_coordinates = np.asarray(x, dtype=np.float64)
    y_coordinates = np.asarray(y, dtype=np.float64)
    if x_coordinates.shape != y_coordinates.shape:
        raise ValueError("x and y coordinate arrays must have matching shapes")
    gravity_x, gravity_y = validate_gravity(gravity)
    return -(gravity_x * x_coordinates + gravity_y * y_coordinates)
