"""Directional approximate Riemann fluxes for the 2D Euler equations."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .eos import sound_speed
from .riemann import RIEMANN_SOLVERS
from .state2d import conservative_to_primitive_2d


def _normal_order(
    conserved: NDArray[np.float64], primitive: NDArray[np.float64], direction: str
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Rotate momentum and velocity components into normal/tangential order."""
    if direction == "x":
        return conserved, primitive
    if direction == "y":
        return conserved[[0, 2, 1, 3]], primitive[[0, 2, 1, 3]]
    raise ValueError("direction must be 'x' or 'y'")


def _physical_flux_normal(
    conserved: NDArray[np.float64], primitive: NDArray[np.float64]
) -> NDArray[np.float64]:
    rho, normal_velocity, tangential_velocity, pressure = primitive
    return np.stack(
        (
            conserved[1],
            rho * normal_velocity**2 + pressure,
            rho * normal_velocity * tangential_velocity,
            normal_velocity * (conserved[3] + pressure),
        ),
        axis=0,
    )


def _restore_order(flux: NDArray[np.float64], direction: str) -> NDArray[np.float64]:
    return flux if direction == "x" else flux[[0, 2, 1, 3]]


def riemann_flux_2d(
    left: ArrayLike,
    right: ArrayLike,
    gamma: float,
    direction: str,
    solver: str = "hll",
) -> NDArray[np.float64]:
    """Return an HLL, HLLC, or Rusanov flux normal to a cell face."""
    if solver not in RIEMANN_SOLVERS:
        raise ValueError(f"unknown Riemann solver {solver!r}; choose from {tuple(RIEMANN_SOLVERS)}")
    left_state = np.asarray(left, dtype=np.float64)
    right_state = np.asarray(right, dtype=np.float64)
    if left_state.shape != right_state.shape or left_state.shape[0] != 4:
        raise ValueError("left and right states must have matching (4, ...) shapes")
    left_primitive = conservative_to_primitive_2d(left_state, gamma)
    right_primitive = conservative_to_primitive_2d(right_state, gamma)
    left_normal, primitive_left_normal = _normal_order(
        left_state, left_primitive, direction
    )
    right_normal, primitive_right_normal = _normal_order(
        right_state, right_primitive, direction
    )
    rho_left, velocity_left, tangent_left, pressure_left = primitive_left_normal
    rho_right, velocity_right, tangent_right, pressure_right = primitive_right_normal
    sound_left = sound_speed(rho_left, pressure_left, gamma)
    sound_right = sound_speed(rho_right, pressure_right, gamma)
    speed_left = np.minimum(velocity_left - sound_left, velocity_right - sound_right)
    speed_right = np.maximum(velocity_left + sound_left, velocity_right + sound_right)
    if np.any(speed_right <= speed_left):
        raise RuntimeError("Riemann wave speeds are not ordered")
    flux_left = _physical_flux_normal(left_normal, primitive_left_normal)
    flux_right = _physical_flux_normal(right_normal, primitive_right_normal)

    if solver == "rusanov":
        maximum_speed = np.maximum(
            np.abs(velocity_left) + sound_left,
            np.abs(velocity_right) + sound_right,
        )
        flux = 0.5 * (flux_left + flux_right) - 0.5 * maximum_speed[np.newaxis, ...] * (
            right_normal - left_normal
        )
        return _restore_order(flux, direction)

    if solver == "hll":
        middle = (
            speed_right[np.newaxis, ...] * flux_left
            - speed_left[np.newaxis, ...] * flux_right
            + (speed_left * speed_right)[np.newaxis, ...]
            * (right_normal - left_normal)
        ) / (speed_right - speed_left)[np.newaxis, ...]
        flux = np.where(
            (speed_left >= 0.0)[np.newaxis, ...],
            flux_left,
            np.where((speed_right <= 0.0)[np.newaxis, ...], flux_right, middle),
        )
        return _restore_order(flux, direction)

    left_mass_speed = rho_left * (speed_left - velocity_left)
    right_mass_speed = rho_right * (speed_right - velocity_right)
    denominator = left_mass_speed - right_mass_speed
    scale = np.maximum(np.abs(left_mass_speed), np.abs(right_mass_speed))
    if np.any(np.abs(denominator) <= np.finfo(float).eps * np.maximum(scale, 1.0)):
        raise RuntimeError("HLLC contact-speed denominator is singular")
    speed_middle = (
        pressure_right
        - pressure_left
        + left_mass_speed * velocity_left
        - right_mass_speed * velocity_right
    ) / denominator
    pressure_star = 0.5 * (
        pressure_left + left_mass_speed * (speed_middle - velocity_left)
        + pressure_right + right_mass_speed * (speed_middle - velocity_right)
    )
    denominator_left = speed_left - speed_middle
    denominator_right = speed_right - speed_middle
    tolerance = np.finfo(float).eps * np.maximum(
        np.maximum(np.abs(speed_left), np.abs(speed_right)), 1.0
    )
    if np.any(np.abs(denominator_left) <= tolerance) or np.any(
        np.abs(denominator_right) <= tolerance
    ):
        raise RuntimeError("HLLC star-state wave speeds are singular")
    rho_star_left = left_mass_speed / denominator_left
    rho_star_right = right_mass_speed / denominator_right
    energy_star_left = (
        (speed_left - velocity_left) * left_normal[3]
        - pressure_left * velocity_left
        + pressure_star * speed_middle
    ) / denominator_left
    energy_star_right = (
        (speed_right - velocity_right) * right_normal[3]
        - pressure_right * velocity_right
        + pressure_star * speed_middle
    ) / denominator_right
    star_left = np.stack(
        (
            rho_star_left,
            rho_star_left * speed_middle,
            rho_star_left * tangent_left,
            energy_star_left,
        ),
        axis=0,
    )
    star_right = np.stack(
        (
            rho_star_right,
            rho_star_right * speed_middle,
            rho_star_right * tangent_right,
            energy_star_right,
        ),
        axis=0,
    )
    flux_star_left = flux_left + speed_left[np.newaxis, ...] * (
        star_left - left_normal
    )
    flux_star_right = flux_right + speed_right[np.newaxis, ...] * (
        star_right - right_normal
    )
    flux = np.where(
        (speed_left >= 0.0)[np.newaxis, ...],
        flux_left,
        np.where(
            (speed_middle >= 0.0)[np.newaxis, ...],
            flux_star_left,
            np.where(
                (speed_right > 0.0)[np.newaxis, ...], flux_star_right, flux_right
            ),
        ),
    )
    return _restore_order(flux, direction)
