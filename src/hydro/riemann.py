"""Approximate Riemann solvers for Euler interface fluxes."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .eos import sound_speed
from .state import conservative_to_primitive, euler_flux


def _validated_states(
    left: ArrayLike, right: ArrayLike, gamma: float
) -> tuple[
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
]:
    """Return conserved/primitive states and Davis signal speeds."""
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
    if np.any(s_right <= s_left):
        raise RuntimeError("Riemann wave speeds are not ordered")
    return u_left, u_right, w_left, w_right, s_left, s_right


def hll_flux(left: ArrayLike, right: ArrayLike, gamma: float) -> NDArray[np.float64]:
    """Compute the HLL flux using the Davis signal-speed estimates."""
    u_left, u_right, _, _, s_left, s_right = _validated_states(
        left, right, gamma
    )
    f_left = euler_flux(u_left, gamma)
    f_right = euler_flux(u_right, gamma)

    denominator = s_right - s_left
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


def rusanov_flux(
    left: ArrayLike, right: ArrayLike, gamma: float
) -> NDArray[np.float64]:
    """Compute the local Lax--Friedrichs (Rusanov) interface flux."""
    u_left, u_right, w_left, w_right, _, _ = _validated_states(left, right, gamma)
    c_left = sound_speed(w_left[0], w_left[2], gamma)
    c_right = sound_speed(w_right[0], w_right[2], gamma)
    maximum_speed = np.maximum(
        np.abs(w_left[1]) + c_left, np.abs(w_right[1]) + c_right
    )
    return 0.5 * (euler_flux(u_left, gamma) + euler_flux(u_right, gamma)) - 0.5 * (
        maximum_speed[np.newaxis, ...] * (u_right - u_left)
    )


def hllc_flux(left: ArrayLike, right: ArrayLike, gamma: float) -> NDArray[np.float64]:
    """Compute the contact-restoring HLLC flux with Davis outer speeds."""
    u_left, u_right, w_left, w_right, s_left, s_right = _validated_states(
        left, right, gamma
    )
    rho_left, velocity_left, pressure_left = w_left
    rho_right, velocity_right, pressure_right = w_right
    f_left = euler_flux(u_left, gamma)
    f_right = euler_flux(u_right, gamma)

    left_mass_speed = rho_left * (s_left - velocity_left)
    right_mass_speed = rho_right * (s_right - velocity_right)
    denominator = left_mass_speed - right_mass_speed
    scale = np.maximum(np.abs(left_mass_speed), np.abs(right_mass_speed))
    if np.any(np.abs(denominator) <= np.finfo(float).eps * np.maximum(scale, 1.0)):
        raise RuntimeError("HLLC contact-speed denominator is singular")
    s_middle = (
        pressure_right
        - pressure_left
        + left_mass_speed * velocity_left
        - right_mass_speed * velocity_right
    ) / denominator

    pressure_star_left = pressure_left + left_mass_speed * (
        s_middle - velocity_left
    )
    pressure_star_right = pressure_right + right_mass_speed * (
        s_middle - velocity_right
    )
    pressure_star = 0.5 * (pressure_star_left + pressure_star_right)

    left_star_denominator = s_left - s_middle
    right_star_denominator = s_right - s_middle
    tolerance = np.finfo(float).eps * np.maximum(
        np.maximum(np.abs(s_left), np.abs(s_right)), 1.0
    )
    if np.any(np.abs(left_star_denominator) <= tolerance) or np.any(
        np.abs(right_star_denominator) <= tolerance
    ):
        raise RuntimeError("HLLC star-state wave speeds are singular")

    rho_star_left = left_mass_speed / left_star_denominator
    rho_star_right = right_mass_speed / right_star_denominator
    energy_star_left = (
        (s_left - velocity_left) * u_left[2]
        - pressure_left * velocity_left
        + pressure_star * s_middle
    ) / left_star_denominator
    energy_star_right = (
        (s_right - velocity_right) * u_right[2]
        - pressure_right * velocity_right
        + pressure_star * s_middle
    ) / right_star_denominator
    u_star_left = np.stack(
        (rho_star_left, rho_star_left * s_middle, energy_star_left), axis=0
    )
    u_star_right = np.stack(
        (rho_star_right, rho_star_right * s_middle, energy_star_right), axis=0
    )
    f_star_left = f_left + s_left[np.newaxis, ...] * (u_star_left - u_left)
    f_star_right = f_right + s_right[np.newaxis, ...] * (u_star_right - u_right)

    return np.where(
        (s_left >= 0.0)[np.newaxis, ...],
        f_left,
        np.where(
            (s_middle >= 0.0)[np.newaxis, ...],
            f_star_left,
            np.where((s_right > 0.0)[np.newaxis, ...], f_star_right, f_right),
        ),
    )


RIEMANN_SOLVERS = {
    "hll": hll_flux,
    "hllc": hllc_flux,
    "rusanov": rusanov_flux,
}


def riemann_flux(
    left: ArrayLike, right: ArrayLike, gamma: float, solver: str = "hll"
) -> NDArray[np.float64]:
    """Dispatch to a configured approximate Riemann flux."""
    try:
        flux_function = RIEMANN_SOLVERS[solver]
    except KeyError as error:
        raise ValueError(
            f"unknown Riemann solver {solver!r}; choose from {tuple(RIEMANN_SOLVERS)}"
        ) from error
    return flux_function(left, right, gamma)
