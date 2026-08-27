import unittest

import numpy as np

from hydro.boundary_conditions import apply_boundaries_2d
from hydro.diagnostics import (
    shock_radius_2d,
    totals_2d,
    transverse_kinetic_energy_2d,
    vorticity_2d,
)
from hydro.exact_riemann import exact_riemann_solution
from hydro.gravity import constant_gravity_source_2d
from hydro.riemann2d import riemann_flux_2d
from hydro.problems import (
    entropy_wave_2d,
    kelvin_helmholtz_2d,
    rayleigh_taylor_2d,
    sedov_taylor_2d,
)
from hydro.reconstruction import reconstruct_interfaces_2d
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

    def test_directional_reconstruction_preserves_uniform_state(self) -> None:
        primitive = np.empty((4, 9, 8))
        primitive[0] = 1.1
        primitive[1] = 0.3
        primitive[2] = -0.2
        primitive[3] = 0.9
        for direction in ("x", "y"):
            with self.subTest(direction=direction):
                left, right = reconstruct_interfaces_2d(primitive, "mc", direction)
                expected_shape = (4, 8, 8) if direction == "x" else (4, 9, 7)
                self.assertEqual(left.shape, expected_shape)
                np.testing.assert_array_equal(left, right)

    def test_periodic_vorticity_matches_smooth_velocity_field(self) -> None:
        resolution = 64
        spacing = 1.0 / resolution
        centers = (np.arange(resolution) + 0.5) * spacing
        x, y = np.meshgrid(centers, centers, indexing="ij")
        primitive = np.empty((4, resolution, resolution))
        primitive[0] = 1.0
        primitive[1] = np.sin(2.0 * np.pi * y)
        primitive[2] = np.cos(2.0 * np.pi * x)
        primitive[3] = 1.0
        numerical = vorticity_2d(primitive, spacing, spacing)
        exact = -2.0 * np.pi * (
            np.sin(2.0 * np.pi * x) + np.cos(2.0 * np.pi * y)
        )
        self.assertLess(float(np.max(np.abs(numerical - exact))), 0.021)

    def test_kelvin_helmholtz_state_and_transverse_energy(self) -> None:
        grid = Grid2D(0.0, 1.0, 64, 0.0, 1.0, 64)
        primitive = kelvin_helmholtz_2d(*grid.mesh)
        self.assertEqual(primitive.shape, (4, 64, 64))
        self.assertGreaterEqual(float(np.min(primitive[0])), 1.0)
        self.assertLessEqual(float(np.max(primitive[0])), 2.0)
        np.testing.assert_allclose(primitive[3], 2.5)
        self.assertGreater(
            transverse_kinetic_energy_2d(primitive, grid.dx * grid.dy), 0.0
        )

    def test_constant_gravity_source_components(self) -> None:
        conserved = np.array([2.0, 0.6, -0.8, 3.0])
        source = constant_gravity_source_2d(conserved, (0.25, -0.5))
        np.testing.assert_allclose(source, [0.0, 0.5, -1.0, 0.55])

    def test_sedov_energy_deposition_is_discretely_normalized(self) -> None:
        grid = Grid2D(-0.5, 0.5, 64, -0.5, 0.5, 64)
        primitive = sedov_taylor_2d(
            *grid.mesh, cell_area=grid.dx * grid.dy, explosion_energy=1.0
        )
        conserved = primitive_to_conservative_2d(primitive, 1.4)
        total_energy = totals_2d(conserved, grid.dx * grid.dy)["energy"]
        ambient_energy = 1.0e-5 / (1.4 - 1.0)
        self.assertAlmostEqual(total_energy - ambient_energy, 1.0, places=13)

    def test_radial_shock_locator_finds_circular_density_jump(self) -> None:
        grid = Grid2D(-0.5, 0.5, 128, -0.5, 0.5, 128)
        x, y = grid.mesh
        radius = np.sqrt(x**2 + y**2)
        primitive = np.empty((4, grid.nx, grid.ny))
        primitive[0] = np.where(radius < 0.2, 4.0, 1.0)
        primitive[1] = 0.0
        primitive[2] = 0.0
        primitive[3] = 1.0
        measured = shock_radius_2d(
            primitive, x, y, grid.dx, minimum_radius=0.1, maximum_radius=0.4
        )
        self.assertLess(abs(measured - 0.2), 2.0 * grid.dx)


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

    def test_muscl_rk2_uniform_periodic_flow_is_preserved(self) -> None:
        grid = Grid2D(0.0, 1.0, 20, 0.0, 1.0, 16)
        solver = Solver2D(
            grid,
            cfl=0.4,
            riemann_solver="hllc",
            reconstruction="muscl",
            limiter="mc",
            integrator="rk2",
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
        solver.run(0.1)
        np.testing.assert_allclose(
            solver.active_conserved, initial, atol=5.0e-15, rtol=0.0
        )

    def test_diagonal_entropy_wave_error_decreases_with_muscl_rk2(self) -> None:
        errors = []
        for resolution in (16, 32):
            grid = Grid2D(0.0, 1.0, resolution, 0.0, 1.0, resolution)
            solver = Solver2D(
                grid,
                cfl=0.4,
                reconstruction="muscl",
                limiter="mc",
                integrator="rk2",
                x_boundary="periodic",
                y_boundary="periodic",
            )
            solver.initialise_function(entropy_wave_2d)
            solver.run(1.0)
            exact = entropy_wave_2d(*grid.mesh, time=1.0)
            errors.append(
                grid.dx * grid.dy * np.sum(np.abs(solver.primitive[0] - exact[0]))
            )
        observed_order = np.log2(errors[0] / errors[1])
        self.assertGreater(observed_order, 1.8)

    def test_muscl_rk2_rotated_sod_is_symmetric_and_positive(self) -> None:
        grid = Grid2D(0.0, 1.0, 40, 0.0, 1.0, 40)
        configuration = dict(
            cfl=0.4,
            riemann_solver="hll",
            reconstruction="muscl",
            limiter="mc",
            integrator="rk2",
        )
        solver_x = Solver2D(grid, **configuration)
        solver_y = Solver2D(grid, **configuration)
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
        self.assertGreater(float(np.min(primitive_x[0])), 0.0)
        self.assertGreater(float(np.min(primitive_x[3])), 0.0)

    def test_uniform_flow_under_gravity_matches_free_fall(self) -> None:
        grid = Grid2D(0.0, 1.0, 12, 0.0, 1.0, 10)
        solver = Solver2D(
            grid,
            cfl=0.4,
            riemann_solver="hllc",
            reconstruction="muscl",
            limiter="mc",
            integrator="rk2",
            gravity=(0.0, -0.5),
            x_boundary="periodic",
            y_boundary="periodic",
        )
        primitive = np.empty((4, grid.nx, grid.ny))
        primitive[0] = 1.2
        primitive[1] = 0.2
        primitive[2] = 0.1
        primitive[3] = 1.0
        solver.initialise(primitive)
        solver.run(0.05)
        final = solver.primitive
        np.testing.assert_allclose(final[0], 1.2, atol=2.0e-14)
        np.testing.assert_allclose(final[1], 0.2, atol=2.0e-14)
        np.testing.assert_allclose(final[2], 0.075, atol=2.0e-14)
        np.testing.assert_allclose(final[3], 1.0, atol=2.0e-14)

    def test_hydrostatic_wall_ghost_pressure_follows_gravity(self) -> None:
        grid = Grid2D(0.0, 1.0, 12, 0.0, 1.0, 12)
        solver = Solver2D(
            grid,
            gravity=(0.0, -0.5),
            x_boundary="periodic",
            y_boundary="hydrostatic",
        )
        solver.initialise_function(
            lambda x, y: rayleigh_taylor_2d(x, y, perturbation_amplitude=0.0)
        )
        full = conservative_to_primitive_2d(solver.conserved, solver.gamma)
        ghost = grid.nghost - 1
        active = grid.nghost
        expected_bottom = -0.5 * (
            full[0, :, ghost] + full[0, :, active]
        ) * solver.gravity[1] * grid.dy
        np.testing.assert_allclose(
            full[3, :, ghost] - full[3, :, active], expected_bottom
        )

