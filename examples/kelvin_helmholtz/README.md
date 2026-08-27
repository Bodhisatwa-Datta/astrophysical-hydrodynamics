# Smooth Kelvin--Helmholtz shear layer

This experiment evolves two smooth, periodic shear layers on the unit square.
The central dense layer moves right, the lighter outer fluid moves left, and a
small sinusoidal transverse velocity perturbation seeds two wavelengths. The
benchmark records density, vorticity, speed, transverse kinetic energy,
conservation, positivity, runtime, and resolution dependence.

From the repository root, I produced the recorded $64^2$, $96^2$, and $128^2$
study with:

```bash
python benchmarks/kelvin_helmholtz.py
```

The run writes the complete time history to
`benchmarks/kelvin_helmholtz_history.csv`, its final summary to
`benchmarks/kelvin_helmholtz_summary.csv`, and both scientific figures under
`figures/`.

All quantities are dimensionless. The calculation is inviscid: changes with
resolution reflect numerical diffusion, not a specified physical viscosity.
