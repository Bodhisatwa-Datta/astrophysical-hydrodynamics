"""Closed-form Sedov--Taylor similarity solution for a uniform 2D medium."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SedovSimilarity2D:
    """Radial cylindrical Sedov solution at one time."""

    radius: np.ndarray
    density: np.ndarray
    radial_velocity: np.ndarray
    pressure: np.ndarray
    shock_radius: float
    shock_speed: float


def _dimensionless_profiles(
    gamma: float,
    samples: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate Kamm's closed-form uniform-density cylindrical solution."""
    if gamma <= 1.0:
        raise ValueError("gamma must exceed one")
    if samples < 128:
        raise ValueError("samples must be at least 128")

    geometry = 2.0
    omega = 0.0
    a = (geometry + 2.0 - omega) * (gamma + 1.0) / 4.0
    b = (gamma + 1.0) / (gamma - 1.0)
    c = (geometry + 2.0 - omega) * gamma / 2.0
    d = (geometry + 2.0 - omega) * (gamma + 1.0) / (
        (geometry + 2.0 - omega) * (gamma + 1.0)
        - 2.0 * (2.0 + geometry * (gamma - 1.0))
    )
    e = (2.0 + geometry * (gamma - 1.0)) / 2.0

    alpha0 = 2.0 / (geometry + 2.0 - omega)
    alpha2 = -(gamma - 1.0) / (
        2.0 * (gamma - 1.0) + geometry - gamma * omega
    )
    alpha1 = (geometry + 2.0 - omega) * gamma / (
        2.0 + geometry * (gamma - 1.0)
    ) * (
        2.0 * (geometry * (2.0 - gamma) - omega)
        / (gamma * (geometry + 2.0 - omega) ** 2)
        - alpha2
    )
    alpha3 = (geometry - omega) / (
        2.0 * (gamma - 1.0) + geometry - gamma * omega
    )
    alpha4 = (
        (geometry + 2.0 - omega)
        * (geometry - omega)
        / (geometry * (2.0 - gamma) - omega)
        * alpha1
    )
    alpha5 = (omega * (1.0 + gamma) - 2.0 * geometry) / (
        geometry * (2.0 - gamma) - omega
    )

    # The standard solution spans V=1/c at the origin to V=1/a at the shock.
    # Geometric spacing in V-1/c resolves the steep density fall near the origin.
    offset = np.geomspace(1.0e-14, 1.0 / a - 1.0 / c, samples)
    similarity_velocity = 1.0 / c + offset
    x1 = a * similarity_velocity
    x2 = b * (c * similarity_velocity - 1.0)
    x3 = d * (1.0 - e * similarity_velocity)
    x4 = b * (1.0 - c * similarity_velocity / gamma)

    radius = x1 ** (-alpha0) * x2 ** (-alpha2) * x3 ** (-alpha1)
    velocity = x1 * radius
    density = (
        x1 ** (alpha0 * omega)
        * x2 ** (alpha3 + alpha2 * omega)
        * x3 ** (alpha4 + alpha1 * omega)
        * x4**alpha5
    )
    pressure = (
        x1 ** (alpha0 * geometry)
        * x3 ** (alpha4 + alpha1 * (omega - 2.0))
        * x4 ** (1.0 + alpha5)
    )

    radius = np.concatenate(([0.0], radius))
    velocity = np.concatenate(([0.0], velocity))
    density = np.concatenate(([0.0], density))
    pressure = np.concatenate(([pressure[0]], pressure))
    return radius, density, velocity, pressure


def sedov_similarity_2d(
    radius: np.ndarray,
    time: float,
    *,
    gamma: float = 1.4,
    explosion_energy: float = 1.0,
    ambient_density: float = 1.0,
    ambient_pressure: float = 0.0,
    samples: int = 8192,
) -> SedovSimilarity2D:
    """Return the strong-shock cylindrical Sedov solution at ``radius``."""
    requested_radius = np.asarray(radius, dtype=float)
    if time <= 0.0:
        raise ValueError("time must be positive")
    if explosion_energy <= 0.0 or ambient_density <= 0.0:
        raise ValueError("explosion energy and ambient density must be positive")
    if np.any(requested_radius < 0.0):
        raise ValueError("radius must be non-negative")

    similarity_radius, density_ratio, velocity_ratio, pressure_ratio = (
        _dimensionless_profiles(gamma, samples)
    )
    postshock_density = (gamma + 1.0) / (gamma - 1.0)
    postshock_velocity = 2.0 / (gamma + 1.0)
    postshock_pressure = 2.0 / (gamma + 1.0)
    dimensionless_density = postshock_density * density_ratio
    dimensionless_velocity = postshock_velocity * velocity_ratio
    dimensionless_pressure = postshock_pressure * pressure_ratio

    energy_density = (
        0.5 * dimensionless_density * dimensionless_velocity**2
        + dimensionless_pressure / (gamma - 1.0)
    )
    radial_energy = energy_density * similarity_radius
    radial_energy_integral = np.sum(
        0.5
        * (radial_energy[:-1] + radial_energy[1:])
        * np.diff(similarity_radius)
    )
    similarity_exponent = 0.5
    normalization = (
        1.0
        / (2.0 * np.pi * similarity_exponent**2 * radial_energy_integral)
    ) ** 0.25
    shock_radius = normalization * (
        explosion_energy / ambient_density
    ) ** 0.25 * np.sqrt(time)
    shock_speed = similarity_exponent * shock_radius / time

    scaled_radius = requested_radius / shock_radius
    inside = scaled_radius <= 1.0
    density = np.full_like(requested_radius, ambient_density)
    radial_velocity = np.zeros_like(requested_radius)
    pressure = np.full_like(requested_radius, ambient_pressure)
    density[inside] = ambient_density * np.interp(
        scaled_radius[inside], similarity_radius, dimensionless_density
    )
    radial_velocity[inside] = shock_speed * np.interp(
        scaled_radius[inside], similarity_radius, dimensionless_velocity
    )
    pressure[inside] = ambient_density * shock_speed**2 * np.interp(
        scaled_radius[inside], similarity_radius, dimensionless_pressure
    )
    return SedovSimilarity2D(
        radius=requested_radius,
        density=density,
        radial_velocity=radial_velocity,
        pressure=pressure,
        shock_radius=float(shock_radius),
        shock_speed=float(shock_speed),
    )
