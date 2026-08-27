# Smooth Rayleigh--Taylor instability

This experiment places a smooth dense layer above a lighter layer under
constant downward gravity. The pressure profile is initialized in analytic
hydrostatic balance. A localized sinusoidal vertical velocity seeds the
instability, while a second calculation with exactly zero perturbation measures
the numerical hydrostatic drift at the same resolution.

I produced the recorded $64^2$, $96^2$, and $128^2$ study with:

```bash
python benchmarks/rayleigh_taylor.py
```

The driver records interface bubble height, spike depth, vertical kinetic
energy, fitted growth rate, vorticity, gas-plus-potential energy, positivity,
runtime, and the unperturbed control history. It writes CSV data under
`benchmarks/` and figures under `figures/`.

All quantities are dimensionless. The current result reaches linear interface
growth at $t=2.5$; it is not presented as a nonlinear mixing-layer study.
