"""Exact ideal-gas Riemann solution for validation of 1D shock tubes.

The implementation follows the pressure-function construction in Toro's
``Riemann Solvers and Numerical Methods for Fluid Dynamics``.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .eos import sound_speed, validate_primitive


def _pressure_function(
    trial: float, rho: float, p: float, sound: float, gamma: float
) -> tuple[float, float]:
    if trial > p:
        a_coeff = 2.0 / ((gamma + 1.0) * rho)
        b_coeff = (gamma - 1.0) / (gamma + 1.0) * p
        root = np.sqrt(a_coeff / (trial + b_coeff))
        value = (trial - p) * root
        derivative = root * (1.0 - 0.5 * (trial - p) / (trial + b_coeff))
    else:
        exponent = (gamma - 1.0) / (2.0 * gamma)
        ratio = trial / p
        value = 2.0 * sound / (gamma - 1.0) * (ratio**exponent - 1.0)
        derivative = ratio ** (-(gamma + 1.0) / (2.0 * gamma)) / (rho * sound)
    return float(value), float(derivative)


def star_region(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
    gamma: float,
    tolerance: float = 1.0e-10,
    max_iterations: int = 100,
) -> tuple[float, float]:
    """Return pressure and velocity in the Riemann problem's star region."""
    rho_l, u_l, p_l = left
    rho_r, u_r, p_r = right
    validate_primitive(rho_l, u_l, p_l)
    validate_primitive(rho_r, u_r, p_r)
    a_l = float(sound_speed(rho_l, p_l, gamma))
    a_r = float(sound_speed(rho_r, p_r, gamma))
    pressure_guess = max(
        tolerance,
        0.5 * (p_l + p_r) - 0.125 * (u_r - u_l) * (rho_l + rho_r) * (a_l + a_r),
    )

    p_star = pressure_guess
    for _ in range(max_iterations):
        f_l, df_l = _pressure_function(p_star, rho_l, p_l, a_l, gamma)
        f_r, df_r = _pressure_function(p_star, rho_r, p_r, a_r, gamma)
        updated = p_star - (f_l + f_r + u_r - u_l) / (df_l + df_r)
        updated = max(tolerance, updated)
        relative_change = 2.0 * abs(updated - p_star) / (updated + p_star)
        p_star = updated
        if relative_change < tolerance:
            break
    else:
        raise RuntimeError("exact Riemann pressure iteration did not converge")

    f_l, _ = _pressure_function(p_star, rho_l, p_l, a_l, gamma)
    f_r, _ = _pressure_function(p_star, rho_r, p_r, a_r, gamma)
    u_star = 0.5 * (u_l + u_r + f_r - f_l)
    return p_star, u_star


def exact_riemann_solution(
    x: ArrayLike,
    time: float,
    left: tuple[float, float, float],
    right: tuple[float, float, float],
    gamma: float,
    discontinuity: float = 0.0,
) -> NDArray[np.float64]:
    """Sample the exact primitive solution ``[rho, u, p]`` at positions ``x``."""
    if time <= 0.0:
        raise ValueError("exact solution requires a positive time")
    coordinates = np.asarray(x, dtype=np.float64)
    p_star, u_star = star_region(left, right, gamma)
    xi_values = (coordinates - discontinuity) / time
    result = np.empty((3,) + coordinates.shape, dtype=np.float64)
    rho_l, u_l, p_l = left
    rho_r, u_r, p_r = right
    a_l = float(sound_speed(rho_l, p_l, gamma))
    a_r = float(sound_speed(rho_r, p_r, gamma))
    g_ratio = (gamma - 1.0) / (gamma + 1.0)

    for index in np.ndindex(coordinates.shape):
        xi = float(xi_values[index])
        if xi <= u_star:
            if p_star > p_l:  # left shock
                speed = u_l - a_l * np.sqrt(
                    (gamma + 1.0) / (2.0 * gamma) * p_star / p_l
                    + (gamma - 1.0) / (2.0 * gamma)
                )
                if xi <= speed:
                    state = left
                else:
                    rho_star = rho_l * (p_star / p_l + g_ratio) / (g_ratio * p_star / p_l + 1.0)
                    state = (rho_star, u_star, p_star)
            else:  # left rarefaction
                head = u_l - a_l
                a_star = a_l * (p_star / p_l) ** ((gamma - 1.0) / (2.0 * gamma))
                tail = u_star - a_star
                if xi <= head:
                    state = left
                elif xi >= tail:
                    state = (rho_l * (p_star / p_l) ** (1.0 / gamma), u_star, p_star)
                else:
                    velocity = 2.0 / (gamma + 1.0) * (a_l + 0.5 * (gamma - 1.0) * u_l + xi)
                    sound = 2.0 / (gamma + 1.0) * (a_l + 0.5 * (gamma - 1.0) * (u_l - xi))
                    ratio = sound / a_l
                    state = (rho_l * ratio ** (2.0 / (gamma - 1.0)), velocity, p_l * ratio ** (2.0 * gamma / (gamma - 1.0)))
        else:
            if p_star > p_r:  # right shock
                speed = u_r + a_r * np.sqrt(
                    (gamma + 1.0) / (2.0 * gamma) * p_star / p_r
                    + (gamma - 1.0) / (2.0 * gamma)
                )
                if xi >= speed:
                    state = right
                else:
                    rho_star = rho_r * (p_star / p_r + g_ratio) / (g_ratio * p_star / p_r + 1.0)
                    state = (rho_star, u_star, p_star)
            else:  # right rarefaction
                head = u_r + a_r
                a_star = a_r * (p_star / p_r) ** ((gamma - 1.0) / (2.0 * gamma))
                tail = u_star + a_star
                if xi >= head:
                    state = right
                elif xi <= tail:
                    state = (rho_r * (p_star / p_r) ** (1.0 / gamma), u_star, p_star)
                else:
                    velocity = 2.0 / (gamma + 1.0) * (-a_r + 0.5 * (gamma - 1.0) * u_r + xi)
                    sound = 2.0 / (gamma + 1.0) * (a_r - 0.5 * (gamma - 1.0) * (u_r - xi))
                    ratio = sound / a_r
                    state = (rho_r * ratio ** (2.0 / (gamma - 1.0)), velocity, p_r * ratio ** (2.0 * gamma / (gamma - 1.0)))
        result[(slice(None),) + index] = state
    return result

