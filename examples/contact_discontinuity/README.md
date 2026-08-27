# Translating contact benchmark

Run from the repository root after installing the package:

```bash
python examples/run_riemann.py contact_discontinuity --cells 400
```

The states are $(\rho,u,p)_L=(1,1,1)$ and
$(\rho,u,p)_R=(0.125,1,1)$ with the contact initially at $x=0.3$. Pressure and
velocity should remain constant while the density jump translates to $x=0.5$
at $t=0.2$. The 1--99% transition width quantifies numerical diffusion.

