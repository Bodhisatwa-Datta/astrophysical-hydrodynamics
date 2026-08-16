# Validation record

## Automated tests

On 2026-08-17 the complete standard-library test suite passed: 11 tests using
Python with NumPy 2.3.5. Coverage includes:

- primitive/conserved round-trip accuracy;
- ideal-gas pressure and sound speed;
- a hand-calculated Euler flux;
- rejection of zero, negative, NaN, and infinite physical states;
- HLL consistency for identical left and right states;
- outflow ghost-cell filling and CFL calculation;
- preservation of a moving uniform state;
- Sod positivity, finite-volume mass/energy conservation, and comparison with
  the internal exact solution.

Reproduce the suite with:

```bash
python -m unittest discover -s tests -v
```

## Sod shock tube

The initial states are

$$
(\rho,u,p)_L=(1,0,1),\qquad
(\rho,u,p)_R=(0.125,0,0.1),\qquad \gamma=1.4,
$$

with the discontinuity at $x=0.5$ in the domain $[0,1]$. The committed figure
was generated with 400 cells, CFL 0.8, and final time 0.2.

| Diagnostic | Measured value |
|---|---:|
| Timesteps | 218 |
| Density $L_1$ error | $6.750936\times10^{-3}$ |
| Relative mass change | $0.0$ (reported precision) |
| Relative total-energy change | $0.0$ (reported precision) |
| Minimum density | $0.125$ |
| Minimum pressure | $0.1$ |
| Maximum Mach number | $0.930420$ |

The numerical solution contains the expected rarefaction, contact, and shock
at positions consistent with the internal exact solution. This statement is
supported by the plotted pointwise comparison and the measured density error;
it is not yet a convergence claim. The first-order HLL method visibly smooths
the rarefaction corners and contact discontinuity.

## Current validation boundary

No experimental convergence order, strong-shock result, contact-advection
metric, rarefaction stress test, or multidimensional benchmark has been run.
Those results must not be inferred from the present Sod test. A multi-resolution
Sod study is the next validation step.
