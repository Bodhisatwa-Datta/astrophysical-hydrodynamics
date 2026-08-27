import math
import unittest

import numpy as np

from hydro.boundary_conditions import apply_periodic
from hydro.problems import SOD
from hydro.reconstruction import minmod, monotonized_central, reconstruct_interfaces, van_leer
from hydro.solver1d import Grid1D, Solver1D
from hydro.validation import run_entropy_wave, run_riemann_problem


class ReconstructionTests(unittest.TestCase):
    def test_limiters_on_monotone_and_opposed_slopes(self) -> None:
        left = np.array([1.0, -1.0, 1.0])
        right = np.array([2.0, -2.0, -1.0])
        np.testing.assert_allclose(minmod(left, right), [1.0, -1.0, 0.0])
        np.testing.assert_allclose(monotonized_central(left, right), [1.5, -1.5, 0.0])
        np.testing.assert_allclose(van_leer(left, right), [4.0 / 3.0, -4.0 / 3.0, 0.0])

    def test_reconstruction_preserves_uniform_primitive_state(self) -> None:
        primitive = np.repeat(np.array([[1.2], [-0.3], [0.7]]), 12, axis=1)
        for limiter in ("minmod", "mc", "vanleer"):
            left, right = reconstruct_interfaces(primitive, limiter)
            np.testing.assert_allclose(left, primitive[:, :-1])
            np.testing.assert_allclose(right, primitive[:, 1:])

    def test_periodic_ghost_cells(self) -> None:
        state = np.zeros((3, 8))
        state[:, 2:6] = np.arange(12).reshape(3, 4)
        apply_periodic(state, 2)
        np.testing.assert_array_equal(state[:, :2], state[:, 4:6])
        np.testing.assert_array_equal(state[:, 6:], state[:, 2:4])


class SecondOrderSolverTests(unittest.TestCase):
    def test_uniform_periodic_state_is_preserved_by_muscl_rk2(self) -> None:
        grid = Grid1D(0.0, 1.0, 48)
        solver = Solver1D(
            grid,
            cfl=0.4,
            reconstruction="muscl",
            limiter="mc",
            integrator="rk2",
            boundary="periodic",
        )
        uniform = np.repeat(np.array([[1.1], [0.4], [0.9]]), grid.n_cells, axis=1)
        solver.initialise(uniform)
        initial = solver.active_conserved.copy()
        solver.run(0.2)
        np.testing.assert_allclose(solver.active_conserved, initial, atol=5.0e-15, rtol=0.0)

    def test_entropy_wave_has_second_order_asymptotic_rate(self) -> None:
        results = [
            run_entropy_wave(n, 0.4, "muscl", "mc", "rk2")
            for n in (32, 64, 128)
        ]
        errors = [float(result.diagnostics["density_L1_error"]) for result in results]
        finest_order = math.log(errors[1] / errors[2], 2.0)
        self.assertGreater(finest_order, 1.8)
        for quantity in ("mass", "momentum", "energy"):
            self.assertLess(
                abs(float(results[-1].diagnostics[f"relative_{quantity}_change"])),
                2.0e-14,
            )

    def test_muscl_sod_solution_remains_positive(self) -> None:
        result = run_riemann_problem(SOD, 120, 0.4, "muscl", "mc", "rk2")
        self.assertGreater(float(result.diagnostics["minimum_density"]), 0.0)
        self.assertGreater(float(result.diagnostics["minimum_pressure"]), 0.0)
