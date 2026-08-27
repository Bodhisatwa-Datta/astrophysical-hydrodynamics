"""Conservative first/second-order unsplit solver for the 2D Euler equations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from .boundary_conditions import apply_boundaries_2d
from .eos import sound_speed, validate_gamma
from .gravity import constant_gravity_source_2d, validate_gravity
from .reconstruction import LIMITERS, reconstruct_interfaces_2d
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
    """Configurable first/second-order unsplit Cartesian Euler solver."""

    def __init__(
        self,
        grid: Grid2D,
        gamma: float = 1.4,
        cfl: float = 0.4,
        riemann_solver: str = "hll",
        x_boundary: str = "outflow",
        y_boundary: str = "outflow",
        reconstruction: str = "constant",
        limiter: str = "mc",
        integrator: str = "euler",
        gravity: tuple[float, float] = (0.0, 0.0),
    ):
        validate_gamma(gamma)
        if not 0.0 < cfl <= 1.0:
            raise ValueError("cfl must lie in (0, 1]")
        if riemann_solver not in RIEMANN_SOLVERS:
            raise ValueError(
                f"unknown Riemann solver {riemann_solver!r}; "
                f"choose from {tuple(RIEMANN_SOLVERS)}"
            )
        valid_x_boundaries = ("outflow", "periodic", "reflective")
        valid_y_boundaries = (*valid_x_boundaries, "hydrostatic")
        if x_boundary not in valid_x_boundaries or y_boundary not in valid_y_boundaries:
            raise ValueError(
                f"x boundary must be chosen from {valid_x_boundaries} and "
                f"y boundary from {valid_y_boundaries}"
            )
        if reconstruction not in ("constant", "muscl"):
            raise ValueError("reconstruction must be 'constant' or 'muscl'")
        if limiter not in LIMITERS:
            raise ValueError(f"unknown limiter {limiter!r}; choose from {tuple(LIMITERS)}")
        if integrator not in ("euler", "rk2"):
            raise ValueError("integrator must be 'euler' or 'rk2'")
        if reconstruction == "muscl" and grid.nghost < 2:
            raise ValueError("MUSCL reconstruction requires at least two ghost cells")
        acceleration = validate_gravity(gravity)
        if y_boundary == "hydrostatic" and acceleration[1] == 0.0:
            raise ValueError("hydrostatic y boundaries require nonzero y gravity")
        self.grid = grid
        self.gamma = gamma
        self.cfl = cfl
        self.riemann_solver = riemann_solver
        self.x_boundary = x_boundary
        self.y_boundary = y_boundary
        self.reconstruction = reconstruction
        self.limiter = limiter
        self.integrator = integrator
        self.gravity = acceleration
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
        hyperbolic_dt = self.cfl / maximum
        acceleration = float(np.hypot(*self.gravity))
        if acceleration == 0.0:
            return hyperbolic_dt
        gravity_dt = np.sqrt(
            self.cfl * min(self.grid.dx, self.grid.dy) / acceleration
        )
        return min(hyperbolic_dt, gravity_dt)

    def step(self, dt: float | None = None) -> float:
        """Advance one Euler or SSP-RK2 unsplit finite-volume step."""
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
        rhs_initial = self._active_rhs(initial)
        if self.integrator == "euler":
            updated = initial[:, x_slice, y_slice] + dt * rhs_initial
        else:
            stage = initial.copy()
            stage[:, x_slice, y_slice] = (
                initial[:, x_slice, y_slice] + dt * rhs_initial
            )
            conservative_to_primitive_2d(stage[:, x_slice, y_slice], self.gamma)
            self._apply_boundaries(stage)
            rhs_stage = self._active_rhs(stage)
            updated = 0.5 * initial[:, x_slice, y_slice] + 0.5 * (
                stage[:, x_slice, y_slice] + dt * rhs_stage
            )
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
        if self.reconstruction == "constant":
            left_x, right_x = state[:, :-1, :], state[:, 1:, :]
            left_y, right_y = state[:, :, :-1], state[:, :, 1:]
        else:
            primitive = conservative_to_primitive_2d(state, self.gamma)
            primitive_left_x, primitive_right_x = reconstruct_interfaces_2d(
                primitive, self.limiter, "x"
            )
            primitive_left_y, primitive_right_y = reconstruct_interfaces_2d(
                primitive, self.limiter, "y"
            )
            left_x = primitive_to_conservative_2d(primitive_left_x, self.gamma)
            right_x = primitive_to_conservative_2d(primitive_right_x, self.gamma)
            left_y = primitive_to_conservative_2d(primitive_left_y, self.gamma)
            right_y = primitive_to_conservative_2d(primitive_right_y, self.gamma)
        flux_x = riemann_flux_2d(
            left_x, right_x, self.gamma, "x", self.riemann_solver
        )
        flux_y = riemann_flux_2d(
            left_y, right_y, self.gamma, "y", self.riemann_solver
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
        rhs = -(divergence_x + divergence_y)
        if self.gravity != (0.0, 0.0):
            rhs += constant_gravity_source_2d(
                state[:, start_x:stop_x, start_y:stop_y], self.gravity
            )
        return rhs

    def _apply_boundaries(self, state: NDArray[np.float64]) -> None:
        y_boundary = "reflective" if self.y_boundary == "hydrostatic" else self.y_boundary
        apply_boundaries_2d(
            state,
            self.grid.nghost,
            self.x_boundary,
            y_boundary,
        )
        if self.y_boundary == "hydrostatic":
            self._apply_hydrostatic_y_pressure(state)

    def _apply_hydrostatic_y_pressure(self, state: NDArray[np.float64]) -> None:
        """Extrapolate ghost pressure with ``dp/dy=rho*g_y`` at solid walls."""
        nghost = self.grid.nghost
        gravity_y = self.gravity[1]
        primitive = conservative_to_primitive_2d(state, self.gamma)
        pressure = primitive[3]
        density = state[0]
        for index in range(nghost - 1, -1, -1):
            pressure[:, index] = pressure[:, index + 1] - 0.5 * (
                density[:, index] + density[:, index + 1]
            ) * gravity_y * self.grid.dy
            self._set_energy_from_pressure(state, pressure[:, index], index)
        top_start = state.shape[2] - nghost
        for index in range(top_start, state.shape[2]):
            pressure[:, index] = pressure[:, index - 1] + 0.5 * (
                density[:, index] + density[:, index - 1]
            ) * gravity_y * self.grid.dy
            self._set_energy_from_pressure(state, pressure[:, index], index)

    def _set_energy_from_pressure(
        self, state: NDArray[np.float64], pressure: NDArray[np.float64], y_index: int
    ) -> None:
        density = state[0, :, y_index]
        kinetic = 0.5 * (
            state[1, :, y_index] ** 2 + state[2, :, y_index] ** 2
        ) / density
        state[3, :, y_index] = pressure / (self.gamma - 1.0) + kinetic

    def _require_initialised(self) -> None:
        if not self._initialised:
            raise RuntimeError("solver has not been initialised")
