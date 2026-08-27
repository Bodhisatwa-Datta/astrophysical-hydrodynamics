"""Conservative first-order unsplit solver for the 2D Euler equations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from .boundary_conditions import apply_boundaries_2d
from .eos import sound_speed, validate_gamma
from .riemann import RIEMANN_SOLVERS
from .riemann2d import riemann_flux_2d
from .state2d import conservative_to_primitive_2d, primitive_to_conservative_2d


@dataclass(frozen=True)
class Grid2D:
    """Uniform Cartesian cell-centred grid with ghost cells."""

    x_min: float
    x_max: float
    nx: int
    y_min: float
    y_max: float
    ny: int
    nghost: int = 2

    def __post_init__(self) -> None:
        bounds = (self.x_min, self.x_max, self.y_min, self.y_max)
        if not all(np.isfinite(value) for value in bounds):
            raise ValueError("grid bounds must be finite")
        if self.x_max <= self.x_min or self.y_max <= self.y_min:
            raise ValueError("upper grid bounds must exceed lower bounds")
        if self.nx < 1 or self.ny < 1 or self.nghost < 1:
            raise ValueError("nx, ny, and nghost must be positive")

    @property
    def dx(self) -> float:
        return (self.x_max - self.x_min) / self.nx

    @property
    def dy(self) -> float:
        return (self.y_max - self.y_min) / self.ny

    @property
    def x_centers(self) -> NDArray[np.float64]:
        return self.x_min + (np.arange(self.nx) + 0.5) * self.dx

    @property
    def y_centers(self) -> NDArray[np.float64]:
        return self.y_min + (np.arange(self.ny) + 0.5) * self.dy

    @property
    def mesh(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        return np.meshgrid(self.x_centers, self.y_centers, indexing="ij")

    @property
    def active(self) -> tuple[slice, slice]:
        return (
            slice(self.nghost, self.nghost + self.nx),
            slice(self.nghost, self.nghost + self.ny),
        )


class Solver2D:
    """First-order unsplit Cartesian finite-volume Euler solver."""

    def __init__(
        self,
        grid: Grid2D,
        gamma: float = 1.4,
        cfl: float = 0.4,
        riemann_solver: str = "hll",
        x_boundary: str = "outflow",
        y_boundary: str = "outflow",
    ):
        validate_gamma(gamma)
        if not 0.0 < cfl <= 1.0:
            raise ValueError("cfl must lie in (0, 1]")
        if riemann_solver not in RIEMANN_SOLVERS:
            raise ValueError(
                f"unknown Riemann solver {riemann_solver!r}; "
                f"choose from {tuple(RIEMANN_SOLVERS)}"
            )
        valid_boundaries = ("outflow", "periodic", "reflective")
        if x_boundary not in valid_boundaries or y_boundary not in valid_boundaries:
            raise ValueError(f"boundaries must be chosen from {valid_boundaries}")
        self.grid = grid
        self.gamma = gamma
        self.cfl = cfl
        self.riemann_solver = riemann_solver
        self.x_boundary = x_boundary
        self.y_boundary = y_boundary
        self.time = 0.0
        self.steps = 0
        self.conserved = np.empty(
            (4, grid.nx + 2 * grid.nghost, grid.ny + 2 * grid.nghost),
            dtype=np.float64,
        )
        self._initialised = False

    @property
    def active_conserved(self) -> NDArray[np.float64]:
        self._require_initialised()
        x_slice, y_slice = self.grid.active
        return self.conserved[:, x_slice, y_slice]

    @property
    def primitive(self) -> NDArray[np.float64]:
        return conservative_to_primitive_2d(self.active_conserved, self.gamma)

    def initialise(self, primitive: NDArray[np.float64]) -> None:
        """Set active cells from an array with shape ``(4, nx, ny)``."""
        values = np.asarray(primitive, dtype=np.float64)
        expected = (4, self.grid.nx, self.grid.ny)
        if values.shape != expected:
            raise ValueError(f"primitive state must have shape {expected}")
        x_slice, y_slice = self.grid.active
        self.conserved[:, x_slice, y_slice] = primitive_to_conservative_2d(
            values, self.gamma
        )
        self._apply_boundaries(self.conserved)
        self.time = 0.0
        self.steps = 0
        self._initialised = True

    def initialise_function(
        self,
        function: Callable[
            [NDArray[np.float64], NDArray[np.float64]], NDArray[np.float64]
        ],
    ) -> None:
        """Initialise from a function evaluated on the active cell-centre mesh."""
        self.initialise(function(*self.grid.mesh))

    def timestep(self) -> float:
        """Return the multidimensional unsplit CFL timestep."""
        primitive = self.primitive
        density, velocity_x, velocity_y, pressure = primitive
        sound = sound_speed(density, pressure, self.gamma)
        inverse_dt = (np.abs(velocity_x) + sound) / self.grid.dx + (
            np.abs(velocity_y) + sound
        ) / self.grid.dy
        maximum = float(np.max(inverse_dt))
        if maximum <= 0.0 or not np.isfinite(maximum):
            raise ValueError("maximum multidimensional signal rate must be positive")
        return self.cfl / maximum

    def step(self, dt: float | None = None) -> float:
        """Advance one first-order unsplit finite-volume step."""
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
        x_slice, y_slice = self.grid.active
        updated = initial[:, x_slice, y_slice] + dt * self._active_rhs(initial)
        conservative_to_primitive_2d(updated, self.gamma)
        self.conserved[:, x_slice, y_slice] = updated
        self._apply_boundaries(self.conserved)
        self.time += dt
        self.steps += 1
        return dt

    def run(self, final_time: float) -> None:
        """Advance to ``final_time``, shortening only the last stable step."""
        self._require_initialised()
        if final_time < self.time or not np.isfinite(final_time):
            raise ValueError("final_time must be finite and no earlier than current time")
        while self.time < final_time:
            self.step(min(self.timestep(), final_time - self.time))

    def _active_rhs(self, state: NDArray[np.float64]) -> NDArray[np.float64]:
        self._apply_boundaries(state)
        flux_x = riemann_flux_2d(
            state[:, :-1, :],
            state[:, 1:, :],
            self.gamma,
            "x",
            self.riemann_solver,
        )
        flux_y = riemann_flux_2d(
            state[:, :, :-1],
            state[:, :, 1:],
            self.gamma,
            "y",
            self.riemann_solver,
        )
        start_x = self.grid.nghost
        stop_x = start_x + self.grid.nx
        start_y = self.grid.nghost
        stop_y = start_y + self.grid.ny
        divergence_x = (
            flux_x[:, start_x:stop_x, start_y:stop_y]
            - flux_x[:, start_x - 1 : stop_x - 1, start_y:stop_y]
        ) / self.grid.dx
        divergence_y = (
            flux_y[:, start_x:stop_x, start_y:stop_y]
            - flux_y[:, start_x:stop_x, start_y - 1 : stop_y - 1]
        ) / self.grid.dy
        return -(divergence_x + divergence_y)

    def _apply_boundaries(self, state: NDArray[np.float64]) -> None:
        apply_boundaries_2d(
            state,
            self.grid.nghost,
            self.x_boundary,
            self.y_boundary,
        )

    def _require_initialised(self) -> None:
        if not self._initialised:
            raise RuntimeError("solver has not been initialised")
