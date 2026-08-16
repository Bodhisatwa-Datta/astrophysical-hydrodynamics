import numpy as np
import unittest

from hydro.diagnostics import totals
from hydro.exact_riemann import exact_riemann_solution
from hydro.problems import sod_initial_condition
from hydro.solver1d import Grid1D, Solver1D


class SolverTests(unittest.TestCase):
    def test_uniform_state_is_preserved(self) -> None:
        grid = Grid1D(0.0, 1.0, 64)
        solver = Solver1D(grid, gamma=1.4, cfl=0.8)
        uniform = np.repeat(np.array([[1.2], [0.3], [0.8]]), grid.n_cells, axis=1)
        solver.initialise(uniform)
        initial = solver.active_conserved.copy()
        solver.run(0.25)
        np.testing.assert_allclose(solver.active_conserved, initial, rtol=0.0, atol=3.0e-15)

    def test_sod_solution_is_positive_conservative_and_close_to_exact(self) -> None:
        grid = Grid1D(0.0, 1.0, 200)
        solver = Solver1D(grid, gamma=1.4, cfl=0.8)
        solver.initialise_function(sod_initial_condition)
        initial_totals = totals(solver.active_conserved, grid.dx)
        solver.run(0.2)
        final_totals = totals(solver.active_conserved, grid.dx)

        primitive = solver.primitive
        self.assertGreater(float(np.min(primitive[0])), 0.0)
        self.assertGreater(float(np.min(primitive[2])), 0.0)
        self.assertLess(abs(final_totals["mass"] - initial_totals["mass"]), 2.0e-14)
        self.assertLess(abs(final_totals["energy"] - initial_totals["energy"]), 2.0e-14)

        exact = exact_riemann_solution(
            grid.centers,
            0.2,
            left=(1.0, 0.0, 1.0),
            right=(0.125, 0.0, 0.1),
            gamma=1.4,
            discontinuity=0.5,
        )
        density_l1 = grid.dx * np.sum(np.abs(primitive[0] - exact[0]))
        self.assertLess(density_l1, 0.03)
