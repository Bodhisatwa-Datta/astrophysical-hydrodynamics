"""TVD slope limiters and primitive-variable MUSCL reconstruction."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .eos import validate_primitive


def _same_shape(*values: ArrayLike) -> tuple[NDArray[np.float64], ...]:
    arrays = tuple(np.asarray(value, dtype=np.float64) for value in values)
    if any(array.shape != arrays[0].shape for array in arrays[1:]):
        raise ValueError("all slope differences must have matching shapes")
    return arrays


def minmod(*values: ArrayLike) -> NDArray[np.float64]:
    """Return the smallest-magnitude argument when all signs agree, else zero."""
    arrays = _same_shape(*values)
    if len(arrays) < 2:
        raise ValueError("minmod requires at least two arguments")
    signs_agree = np.logical_or(
        np.all(np.stack(arrays) > 0.0, axis=0),
        np.all(np.stack(arrays) < 0.0, axis=0),
    )
    magnitude = np.min(np.stack(tuple(np.abs(value) for value in arrays)), axis=0)
    return np.where(signs_agree, np.sign(arrays[0]) * magnitude, 0.0)


def monotonized_central(left: ArrayLike, right: ArrayLike) -> NDArray[np.float64]:
    """Return the monotonized-central (MC) limited slope."""
    left_array, right_array = _same_shape(left, right)
    return minmod(
        0.5 * (left_array + right_array), 2.0 * left_array, 2.0 * right_array
    )


def van_leer(left: ArrayLike, right: ArrayLike) -> NDArray[np.float64]:
    """Return the harmonic van Leer limited slope."""
    left_array, right_array = _same_shape(left, right)
    denominator = left_array + right_array
    with np.errstate(divide="ignore", invalid="ignore"):
        harmonic = 2.0 * left_array * right_array / denominator
    return np.where(left_array * right_array > 0.0, harmonic, 0.0)


LIMITERS: dict[str, Callable[[ArrayLike, ArrayLike], NDArray[np.float64]]] = {
    "minmod": lambda left, right: minmod(left, right),
    "mc": monotonized_central,
    "vanleer": van_leer,
}


def limited_slopes(primitive: ArrayLike, limiter: str) -> NDArray[np.float64]:
    """Compute one limited slope per cell, leaving end-cell slopes at zero."""
    values = np.asarray(primitive, dtype=np.float64)
    if values.ndim != 2 or values.shape[0] != 3:
        raise ValueError("primitive state must have shape (3, n_total)")
    try:
        limiter_function = LIMITERS[limiter]
    except KeyError as error:
        raise ValueError(f"unknown limiter {limiter!r}; choose from {tuple(LIMITERS)}") from error
    slopes = np.zeros_like(values)
    backward = values[:, 1:-1] - values[:, :-2]
    forward = values[:, 2:] - values[:, 1:-1]
    slopes[:, 1:-1] = limiter_function(backward, forward)
    return slopes


def reconstruct_interfaces(
    primitive: ArrayLike, limiter: str
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return left/right primitive states at every cell interface.

    Density, velocity, and pressure are reconstructed directly. TVD limiting
    keeps reconstructed density and pressure between neighboring cell values;
    positivity is nevertheless checked explicitly and never repaired by floors.
    """
    values = np.asarray(primitive, dtype=np.float64)
    if values.ndim != 2 or values.shape[0] != 3:
        raise ValueError("primitive state must have shape (3, n_total)")
    validate_primitive(values[0], values[1], values[2])
    slopes = limited_slopes(values, limiter)
    left = values[:, :-1] + 0.5 * slopes[:, :-1]
    right = values[:, 1:] - 0.5 * slopes[:, 1:]
    validate_primitive(left[0], left[1], left[2])
    validate_primitive(right[0], right[1], right[2])
    return left, right

