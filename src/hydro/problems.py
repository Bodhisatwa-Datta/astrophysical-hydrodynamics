"""Documented initial conditions for hydrodynamic validation problems."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class RiemannProblem:
    """Definition of a two-state ideal-gas Riemann problem."""

    name: str
    left: tuple[float, float, float]
    right: tuple[float, float, float]
    gamma: float
    discontinuity: float
    final_time: float

    def initial_condition(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """Sample the cell-centred primitive initial condition."""
        coordinates = np.asarray(x, dtype=np.float64)
        left_mask = coordinates < self.discontinuity
        return np.stack(
            tuple(
                np.where(left_mask, left_value, right_value)
                for left_value, right_value in zip(self.left, self.right)
            ),
            axis=0,
        )


SOD = RiemannProblem(
    name="sod",
    left=(1.0, 0.0, 1.0),
    right=(0.125, 0.0, 0.1),
    gamma=1.4,
    discontinuity=0.5,
    final_time=0.2,
)

# A deliberately severe pressure-jump test (p_L / p_R = 10^5).  The short
# final time keeps all waves away from the boundaries on [0, 1].
STRONG_SHOCK = RiemannProblem(
    name="strong_shock",
    left=(1.0, 0.0, 1000.0),
    right=(1.0, 0.0, 0.01),
    gamma=1.4,
    discontinuity=0.5,
    final_time=0.01,
)

# Pure contact: pressure and velocity match, while density jumps.  At t=0.2
# the exact discontinuity has translated from x=0.3 to x=0.5.
CONTACT_DISCONTINUITY = RiemannProblem(
    name="contact_discontinuity",
    left=(1.0, 1.0, 1.0),
    right=(0.125, 1.0, 1.0),
    gamma=1.4,
    discontinuity=0.3,
    final_time=0.2,
)

# Symmetric expansion from Toro's standard strong-rarefaction test.  It has a
# very low but non-vacuum central pressure for gamma=1.4.
STRONG_RAREFACTION = RiemannProblem(
    name="rarefaction",
    left=(1.0, -2.0, 0.4),
    right=(1.0, 2.0, 0.4),
    gamma=1.4,
    discontinuity=0.5,
    final_time=0.15,
)

RIEMANN_PROBLEMS = {
    problem.name: problem
    for problem in (SOD, STRONG_SHOCK, CONTACT_DISCONTINUITY, STRONG_RAREFACTION)
}


def entropy_wave(
    x: NDArray[np.float64],
    time: float = 0.0,
    amplitude: float = 0.2,
    velocity: float = 1.0,
    pressure: float = 1.0,
) -> NDArray[np.float64]:
    """Return a smooth periodic entropy wave advected at constant velocity.

    Constant velocity and pressure make this an exact Euler solution. The
    density profile is periodic on a unit domain and returns to its initial
    position after one crossing time when ``velocity=1``.
    """
    coordinates = np.asarray(x, dtype=np.float64)
    phase = 2.0 * np.pi * (coordinates - velocity * time)
    density = 1.0 + amplitude * np.sin(phase)
    return np.stack(
        (
            density,
            np.full_like(coordinates, velocity),
            np.full_like(coordinates, pressure),
        ),
        axis=0,
    )


def entropy_wave_2d(
    x: NDArray[np.float64],
    y: NDArray[np.float64],
    time: float = 0.0,
    amplitude: float = 0.2,
    velocity_x: float = 0.7,
    velocity_y: float = 0.3,
    pressure: float = 1.0,
) -> NDArray[np.float64]:
    """Return an exact diagonal entropy wave on a periodic unit square.

    The density varies along ``x + y`` while pressure and both velocity
    components are constant.  With the default velocities, the profile
    returns to its initial position after ``time=1``.
    """
    x_coordinates = np.asarray(x, dtype=np.float64)
    y_coordinates = np.asarray(y, dtype=np.float64)
    if x_coordinates.shape != y_coordinates.shape:
        raise ValueError("x and y coordinate arrays must have matching shapes")
    phase = 2.0 * np.pi * (
        x_coordinates
        + y_coordinates
        - (velocity_x + velocity_y) * time
    )
    density = 1.0 + amplitude * np.sin(phase)
    return np.stack(
        (
            density,
            np.full_like(density, velocity_x),
            np.full_like(density, velocity_y),
            np.full_like(density, pressure),
        ),
        axis=0,
    )


def kelvin_helmholtz_2d(
    x: NDArray[np.float64],
    y: NDArray[np.float64],
    density_inner: float = 2.0,
    density_outer: float = 1.0,
    shear_velocity: float = 0.5,
    pressure: float = 2.5,
    shear_width: float = 0.025,
    perturbation_amplitude: float = 0.01,
    perturbation_width: float = 0.05,
) -> NDArray[np.float64]:
    """Return a smooth periodic double-shear Kelvin--Helmholtz state.

    The central layer, between y=0.25 and y=0.75, moves to the right while
    the outer fluid moves to the left. A sinusoidal transverse perturbation is
    localised around both shear layers. Coordinates are dimensionless and are
    expected on the unit square.
    """
    x_coordinates = np.asarray(x, dtype=np.float64)
    y_coordinates = np.asarray(y, dtype=np.float64)
    if x_coordinates.shape != y_coordinates.shape:
        raise ValueError("x and y coordinate arrays must have matching shapes")
    positive_parameters = (
        density_inner,
        density_outer,
        pressure,
        shear_width,
        perturbation_width,
    )
    if not all(np.isfinite(value) and value > 0.0 for value in positive_parameters):
        raise ValueError("densities, pressure, and layer widths must be positive")
    if not np.isfinite(shear_velocity) or not np.isfinite(perturbation_amplitude):
        raise ValueError("velocity parameters must be finite")

    layer = 0.5 * (
        np.tanh((y_coordinates - 0.25) / shear_width)
        - np.tanh((y_coordinates - 0.75) / shear_width)
    )
    density = density_outer + (density_inner - density_outer) * layer
    velocity_x = -shear_velocity + 2.0 * shear_velocity * layer

    def periodic_distance(center: float) -> NDArray[np.float64]:
        return (y_coordinates - center + 0.5) % 1.0 - 0.5

    envelope = np.exp(
        -0.5 * (periodic_distance(0.25) / perturbation_width) ** 2
    ) + np.exp(-0.5 * (periodic_distance(0.75) / perturbation_width) ** 2)
    velocity_y = (
        perturbation_amplitude * np.sin(4.0 * np.pi * x_coordinates) * envelope
    )
    return np.stack(
        (
            density,
            velocity_x,
            velocity_y,
            np.full_like(density, pressure),
        ),
        axis=0,
    )


def rayleigh_taylor_2d(
    x: NDArray[np.float64],
    y: NDArray[np.float64],
    gravity_y: float = -0.5,
    density_light: float = 1.0,
    density_heavy: float = 2.0,
    interface: float = 0.5,
    transition_width: float = 0.025,
    interface_pressure: float = 2.5,
    perturbation_amplitude: float = 0.0025,
    perturbation_width: float = 0.05,
) -> NDArray[np.float64]:
    """Return a smooth hydrostatic heavy-over-light Rayleigh--Taylor state."""
    x_coordinates = np.asarray(x, dtype=np.float64)
    y_coordinates = np.asarray(y, dtype=np.float64)
    if x_coordinates.shape != y_coordinates.shape:
        raise ValueError("x and y coordinate arrays must have matching shapes")
    if gravity_y >= 0.0 or not np.isfinite(gravity_y):
        raise ValueError("Rayleigh--Taylor gravity_y must be finite and negative")
    positive_parameters = (
        density_light,
        density_heavy,
        transition_width,
        interface_pressure,
        perturbation_width,
    )
    if not all(np.isfinite(value) and value > 0.0 for value in positive_parameters):
        raise ValueError("densities, pressure, and widths must be positive")
    if density_heavy <= density_light:
        raise ValueError("density_heavy must exceed density_light")
    if not 0.0 < interface < 1.0:
        raise ValueError("interface must lie inside the unit-height domain")

    density_mean = 0.5 * (density_heavy + density_light)
    density_jump = 0.5 * (density_heavy - density_light)
    offset = y_coordinates - interface
    density = density_mean + density_jump * np.tanh(offset / transition_width)
    hydrostatic_integral = (
        density_mean * offset
        + density_jump
        * transition_width
        * np.log(np.cosh(offset / transition_width))
    )
    pressure = interface_pressure + gravity_y * hydrostatic_integral
    velocity_y = (
        perturbation_amplitude
        * np.sin(4.0 * np.pi * x_coordinates)
        * np.exp(-0.5 * (offset / perturbation_width) ** 2)
    )
    return np.stack(
        (density, np.zeros_like(density), velocity_y, pressure), axis=0
    )


def sedov_taylor_2d(
    x: NDArray[np.float64],
    y: NDArray[np.float64],
    cell_area: float,
    gamma: float = 1.4,
    ambient_density: float = 1.0,
    ambient_pressure: float = 1.0e-5,
    explosion_energy: float = 1.0,
    injection_radius: float = 0.05,
    center: tuple[float, float] = (0.0, 0.0),
) -> NDArray[np.float64]:
    """Return a finite-radius, discretely normalized Sedov blast state."""
    x_coordinates = np.asarray(x, dtype=np.float64)
    y_coordinates = np.asarray(y, dtype=np.float64)
    if x_coordinates.shape != y_coordinates.shape:
        raise ValueError("x and y coordinate arrays must have matching shapes")
    positive = (
        cell_area,
        gamma - 1.0,
        ambient_density,
        ambient_pressure,
        explosion_energy,
        injection_radius,
    )
    if not all(np.isfinite(value) and value > 0.0 for value in positive):
        raise ValueError("cell area, thermodynamic values, energy, and radius must be positive")
    center_x, center_y = center
    radius = np.sqrt(
        (x_coordinates - center_x) ** 2 + (y_coordinates - center_y) ** 2
    )
    normalized_radius = radius / injection_radius
    weights = np.where(
        normalized_radius < 1.0,
        (1.0 - normalized_radius**2) ** 2,
        0.0,
    )
    normalization = float(cell_area * np.sum(weights))
    if normalization <= 0.0:
        raise ValueError("injection radius contains no cell centres")
    deposited_energy_density = explosion_energy * weights / normalization
    pressure = ambient_pressure + (gamma - 1.0) * deposited_energy_density
    density = np.full_like(pressure, ambient_density)
    return np.stack(
        (density, np.zeros_like(density), np.zeros_like(density), pressure), axis=0
    )


def sod_initial_condition(
    x: NDArray[np.float64], discontinuity: float = 0.5
) -> NDArray[np.float64]:
    """Return the standard Sod states ``(1,0,1)`` and ``(0.125,0,0.1)``."""
    problem = RiemannProblem(
        name=SOD.name,
        left=SOD.left,
        right=SOD.right,
        gamma=SOD.gamma,
        discontinuity=discontinuity,
        final_time=SOD.final_time,
    )
    return problem.initial_condition(x)
