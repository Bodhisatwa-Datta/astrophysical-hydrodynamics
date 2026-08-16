"""Approximate Riemann solvers for Euler interface fluxes."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .eos import sound_speed
from .state import conservative_to_primitive, euler_flux


def hll_flux(left: ArrayLike, right: ArrayLike, gamma: float) -> NDArray[np.float64]:
    """Compute the HLL flux using the Davis signal-speed estimates."""
    u_left = np.asarray(left, dtype=np.float64)
    u_right = np.asarray(right, dtype=np.float64)
    if u_left.shape != u_right.shape or u_left.shape[0] != 3:
        raise ValueError("left and right states must have matching (3, ...) shapes")

    w_left = conservative_to_primitive(u_left, gamma)
    w_right = conservative_to_primitive(u_right, gamma)
    c_left = sound_speed(w_left[0], w_left[2], gamma)
    c_right = sound_speed(w_right[0], w_right[2], gamma)
    s_left = np.minimum(w_left[1] - c_left, w_right[1] - c_right)
    s_right = np.maximum(w_left[1] + c_left, w_right[1] + c_right)
    f_left = euler_flux(u_left, gamma)
    f_right = euler_flux(u_right, gamma)

    denominator = s_right - s_left
    if np.any(denominator <= 0.0):
        raise RuntimeError("HLL wave speeds are not ordered")
    middle = (
        s_right[np.newaxis, ...] * f_left
        - s_left[np.newaxis, ...] * f_right
        + (s_left * s_right)[np.newaxis, ...] * (u_right - u_left)
    ) / denominator[np.newaxis, ...]
    return np.where(
        (s_left >= 0.0)[np.newaxis, ...],
        f_left,
        np.where((s_right <= 0.0)[np.newaxis, ...], f_right, middle),
    )

