"""Documented initial conditions for one-dimensional validation problems."""

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
