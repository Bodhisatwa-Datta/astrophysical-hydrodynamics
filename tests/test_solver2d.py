import unittest

import numpy as np

from hydro.boundary_conditions import apply_boundaries_2d
from hydro.diagnostics import totals_2d
from hydro.exact_riemann import exact_riemann_solution
from hydro.riemann2d import riemann_flux_2d
from hydro.solver2d import Grid2D, Solver2D
from hydro.state2d import (
    conservative_to_primitive_2d,
    euler_flux_2d,
    primitive_to_conservative_2d,
)


def sod_x(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    left = x < 0.5
    return np.stack(
        (
            np.where(left, 1.0, 0.125),
            np.zeros_like(x),
            np.zeros_like(x),
            np.where(left, 1.0, 0.1),
        ),
        axis=0,
    )


def sod_y(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    lower = y < 0.5
    return np.stack(
        (
            np.where(lower, 1.0, 0.125),
            np.zeros_like(y),
            np.zeros_like(y),
            np.where(lower, 1.0, 0.1),
        ),
        axis=0,
    )


class State2DTests(unittest.TestCase):
    def test_primitive_conservative_round_trip(self) -> None:
        primitive = np.array(
            [
                [[1.0, 0.7], [1.2, 0.9]],
                [[0.2, -0.3], [0.1, 0.0]],
                [[-0.1, 0.4], [0.2, -0.2]],
                [[1.0, 0.6], [0.8, 1.1]],
            ]
        )
        conserved = primitive_to_conservative_2d(primitive, 1.4)
        recovered = conservative_to_primitive_2d(conserved, 1.4)
        np.testing.assert_allclose(recovered, primitive, atol=2.0e-15, rtol=2.0e-15)

    def test_all_directional_fluxes_are_consistent(self) -> None:
        primitive = np.array([1.0, 0.3, -0.2, 1.0])
        conserved = primitive_to_conservative_2d(primitive, 1.4)
        for direction in ("x", "y"):
            physical = euler_flux_2d(conserved, 1.4, direction)
            for solver in ("hll", "hllc", "rusanov"):
                with self.subTest(direction=direction, solver=solver):
                    numerical = riemann_flux_2d(
                        conserved, conserved, 1.4, direction, solver
                    )
                    np.testing.assert_allclose(numerical, physical, atol=2.0e-15)


class Boundary2DTests(unittest.TestCase):
    def test_reflective_boundaries_reverse_only_normal_momentum(self) -> None:
        state = np.zeros((4, 8, 7))
        state[:, 2:6, 2:5] = np.arange(48).reshape(4, 4, 3) + 1.0
        apply_boundaries_2d(state, 2, "reflective", "reflective")

        np.testing.assert_array_equal(state[0, :2, 2:5], state[0, 2:4, 2:5][::-1])
        np.testing.assert_array_equal(state[1, :2, 2:5], -state[1, 2:4, 2:5][::-1])
        np.testing.assert_array_equal(state[2, :2, 2:5], state[2, 2:4, 2:5][::-1])
        np.testing.assert_array_equal(state[0, 2:6, :2], state[0, 2:6, 2:4][:, ::-1])
        np.testing.assert_array_equal(state[2, 2:6, :2], -state[2, 2:6, 2:4][:, ::-1])


class Solver2DTests(unittest.TestCase):
    def test_uniform_periodic_flow_is_preserved_and_conservative(self) -> None:
        grid = Grid2D(0.0, 1.0, 20, 0.0, 1.0, 16)
        for flux in ("hll", "hllc", "rusanov"):
            with self.subTest(flux=flux):
                solver = Solver2D(
                    grid,
                    cfl=0.4,
                    riemann_solver=flux,
                    x_boundary="periodic",
                    y_boundary="periodic",
                )
                uniform = np.empty((4, grid.nx, grid.ny))
                uniform[0] = 1.1
                uniform[1] = 0.3
                uniform[2] = -0.2
                uniform[3] = 0.9
                solver.initialise(uniform)
                initial = solver.active_conserved.copy()
                initial_totals = totals_2d(initial, grid.dx * grid.dy)
                solver.run(0.1)
                np.testing.assert_allclose(
                    solver.active_conserved, initial, atol=5.0e-15, rtol=0.0
                )
                final_totals = totals_2d(
                    solver.active_conserved, grid.dx * grid.dy
                )
                for quantity in initial_totals:
                    self.assertAlmostEqual(
                        final_totals[quantity], initial_totals[quantity], places=14
                    )

    def test_multidimensional_cfl_rate(self) -> None:
        grid = Grid2D(0.0, 1.0, 10, 0.0, 1.0, 20)
        solver = Solver2D(
            grid, cfl=0.4, x_boundary="periodic", y_boundary="periodic"
        )
        primitive = np.empty((4, grid.nx, grid.ny))
        primitive[0] = 1.0
        primitive[1] = 0.3
        primitive[2] = -0.2
        primitive[3] = 1.0
        solver.initialise(primitive)
        sound = np.sqrt(1.4)
        expected = 0.4 / ((0.3 + sound) / grid.dx + (0.2 + sound) / grid.dy)
        self.assertAlmostEqual(solver.timestep(), expected)

    def test_rotated_sod_solutions_match(self) -> None:
        grid = Grid2D(0.0, 1.0, 48, 0.0, 1.0, 48)
        solver_x = Solver2D(grid, cfl=0.4, riemann_solver="hll")
        solver_y = Solver2D(grid, cfl=0.4, riemann_solver="hll")
        solver_x.initialise_function(sod_x)
        solver_y.initialise_function(sod_y)
        solver_x.run(0.08)
        solver_y.run(0.08)
        primitive_x = solver_x.primitive
        primitive_y = solver_y.primitive
        np.testing.assert_allclose(primitive_x[0], primitive_y[0].T, atol=3.0e-15)
        np.testing.assert_allclose(primitive_x[1], primitive_y[2].T, atol=3.0e-15)
        np.testing.assert_allclose(primitive_x[2], primitive_y[1].T, atol=3.0e-15)
        np.testing.assert_allclose(primitive_x[3], primitive_y[3].T, atol=3.0e-15)

        exact = exact_riemann_solution(
            grid.x_centers,
            0.08,
            left=(1.0, 0.0, 1.0),
            right=(0.125, 0.0, 0.1),
            gamma=1.4,
            discontinuity=0.5,
        )
        density_slice = primitive_x[0, :, grid.ny // 2]
        density_l1 = grid.dx * np.sum(np.abs(density_slice - exact[0]))
        self.assertLess(density_l1, 0.03)

