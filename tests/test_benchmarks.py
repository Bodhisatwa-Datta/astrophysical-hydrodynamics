import unittest

import numpy as np

from hydro.problems import CONTACT_DISCONTINUITY, SOD, STRONG_RAREFACTION, STRONG_SHOCK
from hydro.validation import contact_transition_width, run_riemann_problem


class RiemannBenchmarkTests(unittest.TestCase):
    def test_sod_density_error_decreases_under_refinement(self) -> None:
        errors = [
            float(run_riemann_problem(SOD, cells).diagnostics["density_L1_error"])
            for cells in (50, 100, 200)
        ]
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])

    def test_contact_preserves_pressure_and_velocity(self) -> None:
        result = run_riemann_problem(CONTACT_DISCONTINUITY, 160)
        np.testing.assert_allclose(result.numerical[1], 1.0, atol=2.0e-14, rtol=0.0)
        np.testing.assert_allclose(result.numerical[2], 1.0, atol=3.0e-14, rtol=0.0)
        self.assertLess(
            abs(float(result.diagnostics["mass_boundary_budget_residual"])), 2.0e-14
        )
        self.assertLess(
            abs(float(result.diagnostics["energy_boundary_budget_residual"])), 2.0e-14
        )
        width = contact_transition_width(result)
        self.assertGreater(width, 0.0)
        # First-order HLL is intentionally diffusive at a contact; this broad
        # bound catches pathological domain-wide smearing without pretending
        # the baseline method resolves the discontinuity sharply.
        self.assertLess(width, 0.25)

    def test_strong_shock_remains_positive_and_matches_exact_structure(self) -> None:
        result = run_riemann_problem(STRONG_SHOCK, 160, cfl=0.7)
        self.assertGreater(float(result.diagnostics["minimum_density"]), 0.0)
        self.assertGreater(float(result.diagnostics["minimum_pressure"]), 0.0)
        self.assertLess(float(result.diagnostics["density_L1_error"]), 0.2)

    def test_strong_rarefaction_remains_positive(self) -> None:
        result = run_riemann_problem(STRONG_RAREFACTION, 160, cfl=0.7)
        self.assertGreater(float(result.diagnostics["minimum_density"]), 0.0)
        self.assertGreater(float(result.diagnostics["minimum_pressure"]), 0.0)
        self.assertLess(float(result.diagnostics["density_L1_error"]), 0.1)
