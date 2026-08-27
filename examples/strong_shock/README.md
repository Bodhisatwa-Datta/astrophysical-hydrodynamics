# Strong pressure-jump benchmark

Run from the repository root after installing the package:

```bash
python examples/run_riemann.py strong_shock --cells 400 --cfl 0.7
```

The states are $(\rho,u,p)_L=(1,0,1000)$ and
$(\rho,u,p)_R=(1,0,0.01)$ with $\gamma=1.4$. The run ends at $t=0.01$ before
the physical waves reach the outflow boundaries. The driver prints exact-error,
positivity, runtime, and boundary-budget diagnostics and writes the figure to
`figures/strong_shock_hll_first_order.png`.

