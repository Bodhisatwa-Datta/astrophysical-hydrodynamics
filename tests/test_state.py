import numpy as np
import unittest

from hydro.eos import UnphysicalStateError, pressure, sound_speed
from hydro.state import conservative_to_primitive, euler_flux, primitive_to_conservative


class StateTests(unittest.TestCase):
    def test_primitive_conservative_round_trip(self) -> None:
        primitive = np.array([[1.0, 0.4, 2.0], [-0.5, 0.0, 1.25], [0.8, 2.5, 0.3]])
        conserved = primitive_to_conservative(primitive, gamma=1.4)
        recovered = conservative_to_primitive(conserved, gamma=1.4)
        np.testing.assert_allclose(recovered, primitive, rtol=2.0e-15, atol=2.0e-15)

    def test_pressure_and_sound_speed(self) -> None:
        primitive = np.array([2.0, 3.0, 5.0])
        conserved = primitive_to_conservative(primitive, gamma=1.4)
        self.assertAlmostEqual(float(pressure(conserved, gamma=1.4)), 5.0)
        self.assertAlmostEqual(float(sound_speed(2.0, 5.0, gamma=1.4)), np.sqrt(3.5))

    def test_physical_flux(self) -> None:
        conserved = primitive_to_conservative(np.array([1.0, 2.0, 1.0]), gamma=1.4)
        np.testing.assert_allclose(euler_flux(conserved, gamma=1.4), [2.0, 5.0, 11.0])

    def test_unphysical_primitive_states_raise(self) -> None:
        invalid_states = (
            np.array([0.0, 0.0, 1.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([1.0, np.nan, 1.0]),
            np.array([1.0, 0.0, np.inf]),
        )
        for primitive in invalid_states:
            with self.subTest(primitive=primitive), self.assertRaises(UnphysicalStateError):
                primitive_to_conservative(primitive, gamma=1.4)

    def test_negative_internal_energy_raises(self) -> None:
        with self.assertRaises(UnphysicalStateError):
            conservative_to_primitive(np.array([1.0, 10.0, 1.0]), gamma=1.4)

    def test_invalid_gamma_raises(self) -> None:
        with self.assertRaises(ValueError):
            primitive_to_conservative(np.array([1.0, 0.0, 1.0]), gamma=1.0)
