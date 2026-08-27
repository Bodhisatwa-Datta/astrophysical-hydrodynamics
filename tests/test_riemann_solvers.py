import unittest

from hydro.problems import CONTACT_DISCONTINUITY, STRONG_SHOCK
from hydro.validation import contact_transition_width, run_riemann_problem


class RiemannSolverComparisonTests(unittest.TestCase):
    def test_all_fluxes_remain_positive_on_strong_shock(self) -> None:
        for flux in ("hll", "hllc", "rusanov"):
            with self.subTest(flux=flux):
                result = run_riemann_problem(
                    STRONG_SHOCK, 120, 0.4, "muscl", "mc", "rk2", flux
                )
                self.assertGreater(float(result.diagnostics["minimum_density"]), 0.0)
                self.assertGreater(float(result.diagnostics["minimum_pressure"]), 0.0)

    def test_hllc_resolves_translating_contact_more_sharply(self) -> None:
        widths = {}
        for flux in ("hll", "hllc", "rusanov"):
            result = run_riemann_problem(
                CONTACT_DISCONTINUITY,
                160,
                0.4,
                "constant",
                "mc",
                "euler",
                flux,
            )
            widths[flux] = contact_transition_width(result)
        self.assertLess(widths["hllc"], widths["hll"])
        self.assertLess(widths["hll"], widths["rusanov"])

