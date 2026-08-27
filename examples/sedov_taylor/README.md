# Two-dimensional Sedov--Taylor blast

This experiment deposits exactly one dimensionless unit of thermal energy in a
compact region of radius 0.05 inside a uniform medium. The deposition kernel is
renormalized on every grid, avoiding a resolution-dependent explosion energy
and avoiding a singular single-cell injection.

I produced the recorded $64^2$, $96^2$, and $128^2$ study with:

```bash
python benchmarks/sedov_taylor.py
```

The driver records shock radius versus time, radial density, velocity, and
pressure profiles, conservation, positivity, runtime, and angular shock-radius
scatter. It fits $R=Ct^\alpha$ and compares the trajectory and full radial
structure with the closed-form cylindrical Sedov solution. For these
parameters the exact shock radius at $t=0.05$ is 0.22451; the measured error
decreases from 4.86% at $64^2$ to 2.26% at $128^2$.

The CSV data are written under `benchmarks/` and the field and radial-analysis
figures under `figures/`. All quantities are dimensionless.
