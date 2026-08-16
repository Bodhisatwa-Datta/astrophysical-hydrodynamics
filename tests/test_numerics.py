import numpy as np
import unittest

from hydro.boundary_conditions import apply_outflow
from hydro.riemann import hll_flux
from hydro.state import euler_flux, primitive_to_conservative
from hydro.timestepping import cfl_timestep


class NumericsTests(unittest.TestCase):
    def test_hll_consistency_for_identical_states(self) -> None:
        primitive = np.array([[1.0, 0.7], [0.2, -0.1], [1.0, 0.4]])
        conserved = primitive_to_conservative(primitive, gamma=1.4)
        np.testing.assert_allclose(hll_flux(conserved, conserved, 1.4), euler_flux(conserved, 1.4))

    def test_outflow_boundary_copies_edge_cells(self) -> None:
        state = np.zeros((3, 8))
        state[:, 2:6] = np.arange(12).reshape(3, 4)
        apply_outflow(state, nghost=2)
        np.testing.assert_array_equal(state[:, :2], np.repeat(state[:, 2:3], 2, axis=1))
        np.testing.assert_array_equal(state[:, -2:], np.repeat(state[:, 5:6], 2, axis=1))

    def test_cfl_timestep_for_stationary_uniform_state(self) -> None:
        gamma = 1.4
        primitive = np.repeat(np.array([[1.0], [0.0], [1.0]]), 10, axis=1)
        conserved = primitive_to_conservative(primitive, gamma)
        expected = 0.8 * 0.1 / np.sqrt(gamma)
        self.assertAlmostEqual(cfl_timestep(conserved, 0.1, gamma, 0.8), expected)
