"""First-order finite-volume solver for the one-dimensional Euler equations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from .boundary_conditions import apply_outflow, apply_periodic
from .reconstruction import LIMITERS, reconstruct_interfaces
from .riemann import RIEMANN_SOLVERS, riemann_flux
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
    """Configurable first/second-order finite-volume Euler solver."""

    def __init__(
        self,
        grid: Grid1D,
        gamma: float = 1.4,
        cfl: float = 0.8,
        reconstruction: str = "constant",
        limiter: str = "mc",
        integrator: str = "euler",
        boundary: str = "outflow",
        riemann_solver: str = "hll",
    ):
        if gamma <= 1.0 or not np.isfinite(gamma):
            raise ValueError("gamma must be finite and greater than one")
        if not 0.0 < cfl <= 1.0:
            raise ValueError("cfl must lie in (0, 1]")
        if reconstruction not in ("constant", "muscl"):
            raise ValueError("reconstruction must be 'constant' or 'muscl'")
        if limiter not in LIMITERS:
            raise ValueError(f"unknown limiter {limiter!r}; choose from {tuple(LIMITERS)}")
        if integrator not in ("euler", "rk2"):
            raise ValueError("integrator must be 'euler' or 'rk2'")
        if boundary not in ("outflow", "periodic"):
            raise ValueError("boundary must be 'outflow' or 'periodic'")
        if riemann_solver not in RIEMANN_SOLVERS:
            raise ValueError(
                f"unknown Riemann solver {riemann_solver!r}; "
                f"choose from {tuple(RIEMANN_SOLVERS)}"
            )
        if reconstruction == "muscl" and grid.nghost < 2:
            raise ValueError("MUSCL reconstruction requires at least two ghost cells")
        self.grid = grid
        self.gamma = gamma
        self.cfl = cfl
        self.reconstruction = reconstruction
        self.limiter = limiter
        self.integrator = integrator
        self.boundary = boundary
        self.riemann_solver = riemann_solver
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
        self._apply_boundary(self.conserved)
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

        initial = self.conserved.copy()
        rhs_initial = self._active_rhs(initial)
        start = self.grid.nghost
        stop = start + self.grid.n_cells

        if self.integrator == "euler":
            updated = initial[:, start:stop] + dt * rhs_initial
        else:
            stage = initial.copy()
            stage[:, start:stop] = initial[:, start:stop] + dt * rhs_initial
            conservative_to_primitive(stage[:, start:stop], self.gamma)
            self._apply_boundary(stage)
            rhs_stage = self._active_rhs(stage)
            updated = 0.5 * initial[:, start:stop] + 0.5 * (
                stage[:, start:stop] + dt * rhs_stage
            )

        self.conserved[:, start:stop] = updated

        # Validate immediately; no density or pressure floors are applied.
        conservative_to_primitive(self.conserved[:, start:stop], self.gamma)
        self._apply_boundary(self.conserved)
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

    def _apply_boundary(self, state: NDArray[np.float64]) -> None:
        if self.boundary == "outflow":
            apply_outflow(state, self.grid.nghost)
        else:
            apply_periodic(state, self.grid.nghost)

    def _active_rhs(self, state: NDArray[np.float64]) -> NDArray[np.float64]:
        """Return the semi-discrete flux divergence in active cells."""
        self._apply_boundary(state)
        if self.reconstruction == "constant":
            interface_flux = riemann_flux(
                state[:, :-1],
                state[:, 1:],
                self.gamma,
                self.riemann_solver,
            )
        else:
            primitive = conservative_to_primitive(state, self.gamma)
            left_primitive, right_primitive = reconstruct_interfaces(
                primitive, self.limiter
            )
            left = primitive_to_conservative(left_primitive, self.gamma)
            right = primitive_to_conservative(right_primitive, self.gamma)
            interface_flux = riemann_flux(
                left, right, self.gamma, self.riemann_solver
            )

        start = self.grid.nghost
        stop = start + self.grid.n_cells
        flux_right = interface_flux[:, start:stop]
        flux_left = interface_flux[:, start - 1 : stop - 1]
        return -(flux_right - flux_left) / self.grid.dx
