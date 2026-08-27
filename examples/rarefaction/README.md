# Strong rarefaction benchmark

Run from the repository root after installing the package:

```bash
python examples/run_riemann.py rarefaction --cells 400 --cfl 0.7
```

The symmetric states are $(\rho,u,p)_L=(1,-2,0.4)$ and
$(\rho,u,p)_R=(1,2,0.4)$ with $\gamma=1.4$, evolved to $t=0.15$. This test
probes positivity in the low-density, low-pressure centre without applying
thermodynamic floors.
