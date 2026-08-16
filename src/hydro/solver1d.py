"""First-order finite-volume solver for the one-dimensional Euler equations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from .boundary_conditions import apply_outflow
from .riemann import hll_flux
from .state import conservative_to_primitive, primitive_to_conservative
from .timestepping import cfl_timestep


@dataclass(frozen=True)
class Grid1D:
    """Uniform cell-centred grid with ghost cells."""

    x_min: float
    x_max: float
    n_cells: int
    nghost: int = 2

    def __post_init__(self) -> None:
        if not np.isfinite(self.x_min) or not np.isfinite(self.x_max):
            raise ValueError("grid bounds must be finite")
        if self.x_max <= self.x_min:
            raise ValueError("x_max must exceed x_min")
        if self.n_cells < 1 or self.nghost < 1:
            raise ValueError("n_cells and nghost must be positive")

    @property
    def dx(self) -> float:
        return (self.x_max - self.x_min) / self.n_cells

    @property
    def centers(self) -> NDArray[np.float64]:
        return self.x_min + (np.arange(self.n_cells) + 0.5) * self.dx

    @property
    def active(self) -> slice:
        return slice(self.nghost, self.nghost + self.n_cells)


class Solver1D:
    """Godunov finite-volume method with piecewise-constant states and HLL fluxes."""

    def __init__(self, grid: Grid1D, gamma: float = 1.4, cfl: float = 0.8):
        if gamma <= 1.0 or not np.isfinite(gamma):
            raise ValueError("gamma must be finite and greater than one")
        if not 0.0 < cfl <= 1.0:
            raise ValueError("cfl must lie in (0, 1]")
        self.grid = grid
        self.gamma = gamma
        self.cfl = cfl
        self.time = 0.0
        self.steps = 0
        self.conserved = np.empty((3, grid.n_cells + 2 * grid.nghost))
        self._initialised = False

    @property
    def active_conserved(self) -> NDArray[np.float64]:
        self._require_initialised()
        return self.conserved[:, self.grid.active]

    @property
    def primitive(self) -> NDArray[np.float64]:
        return conservative_to_primitive(self.active_conserved, self.gamma)

    def initialise(self, primitive: NDArray[np.float64]) -> None:
        """Set active cells from an array with shape ``(3, n_cells)``."""
        values = np.asarray(primitive, dtype=np.float64)
        if values.shape != (3, self.grid.n_cells):
            raise ValueError(f"primitive state must have shape (3, {self.grid.n_cells})")
        self.conserved[:, self.grid.active] = primitive_to_conservative(values, self.gamma)
        apply_outflow(self.conserved, self.grid.nghost)
        self.time = 0.0
        self.steps = 0
        self._initialised = True

    def initialise_function(
        self, function: Callable[[NDArray[np.float64]], NDArray[np.float64]]
    ) -> None:
        """Initialise from a function evaluated at active cell centres."""
        self.initialise(function(self.grid.centers))

    def timestep(self) -> float:
        self._require_initialised()
        return cfl_timestep(
            self.conserved, self.grid.dx, self.gamma, self.cfl, self.grid.nghost
        )

    def step(self, dt: float | None = None) -> float:
        """Advance one explicit Euler step and return the timestep used."""
        self._require_initialised()
        stable_dt = self.timestep()
        if dt is None:
            dt = stable_dt
        if dt <= 0.0 or not np.isfinite(dt):
            raise ValueError("dt must be finite and positive")
        tolerance = 16.0 * np.finfo(float).eps * stable_dt
        if dt > stable_dt + tolerance:
            raise ValueError(f"requested dt={dt:g} exceeds CFL limit {stable_dt:g}")

        apply_outflow(self.conserved, self.grid.nghost)
        interface_flux = hll_flux(
            self.conserved[:, :-1], self.conserved[:, 1:], self.gamma
        )
        start = self.grid.nghost
        stop = start + self.grid.n_cells
        flux_right = interface_flux[:, start:stop]
        flux_left = interface_flux[:, start - 1 : stop - 1]
        self.conserved[:, start:stop] -= dt / self.grid.dx * (flux_right - flux_left)

        # Validate immediately; no density or pressure floors are applied.
        conservative_to_primitive(self.conserved[:, start:stop], self.gamma)
        apply_outflow(self.conserved, self.grid.nghost)
        self.time += dt
        self.steps += 1
        return dt

    def run(self, final_time: float) -> None:
        """Advance to ``final_time``, shortening only the final stable step."""
        self._require_initialised()
        if final_time < self.time or not np.isfinite(final_time):
            raise ValueError("final_time must be finite and no earlier than current time")
        while self.time < final_time:
            dt = min(self.timestep(), final_time - self.time)
            self.step(dt)

    def _require_initialised(self) -> None:
        if not self._initialised:
            raise RuntimeError("solver has not been initialised")

