# Validation record

## Validation procedure

I validated the solver in stages rather than treating a successful shock-tube
plot as sufficient evidence. I began with state conversion, the equation of
state, sound speed, and a hand-calculated physical flux. I then checked that
the HLL flux reduces to the physical Euler flux when both interface states are
identical, that the outflow ghost cells contain the intended extrapolated
state, and that a spatially uniform moving fluid remains unchanged.

After those local checks, I evolved the Sod, strong-pressure-jump, translating
contact, and symmetric-rarefaction problems. For each run I inspected density,
velocity, pressure, and specific internal energy against the internally
computed exact Riemann solution. I also recorded minimum density and pressure,
integrated errors, and conservation residuals after accounting for fluxes
through the open boundaries.

The small regression checks used during development can be repeated with:

```bash
python -m unittest discover -s tests -v
```

## Sod shock tube

The initial states are

$$
(\rho,u,p)_L=(1,0,1),\qquad
(\rho,u,p)_R=(0.125,0,0.1),\qquad \gamma=1.4,
$$

with the discontinuity at $x=0.5$ in $[0,1]$. I ran the reference case with
400 cells, CFL 0.8, and final time 0.2, then overlaid the numerical and exact
profiles before calculating the errors below.

| Diagnostic | Measured value |
|---|---:|
| Timesteps | 218 |
| Density $L_1$ error | $6.750936\times10^{-3}$ |
| Relative mass change | $0.0$ (reported precision) |
| Relative total-energy change | $0.0$ (reported precision) |
| Minimum density | $0.125$ |
| Minimum pressure | $0.1$ |

The rarefaction, contact, and shock positions agree with the internal exact
solution at the plotted resolution. The first-order HLL method visibly smooths
the rarefaction corners and contact discontinuity.

## Sod resolution study

I repeated the same Sod calculation at 50, 100, 200, 400, and 800 cells without
changing the CFL number or final time. Complete machine-readable results are in
`benchmarks/convergence/sod_first_order.csv`.

| Cells | Density $L_1$ | Order | Velocity $L_1$ | Order | Pressure $L_1$ | Order |
|---:|---:|---:|---:|---:|---:|---:|
| 50  | $2.485506\times10^{-2}$ | --    | $4.212150\times10^{-2}$ | --    | $2.214574\times10^{-2}$ | --    |
| 100 | $1.675971\times10^{-2}$ | 0.569 | $2.545437\times10^{-2}$ | 0.727 | $1.350572\times10^{-2}$ | 0.713 |
| 200 | $1.065236\times10^{-2}$ | 0.654 | $1.459752\times10^{-2}$ | 0.802 | $8.105195\times10^{-3}$ | 0.737 |
| 400 | $6.750936\times10^{-3}$ | 0.658 | $8.298003\times10^{-3}$ | 0.815 | $4.789993\times10^{-3}$ | 0.759 |
| 800 | $4.270267\times10^{-3}$ | 0.661 | $4.752661\times10^{-3}$ | 0.804 | $2.805447\times10^{-3}$ | 0.772 |

All errors decrease monotonically. The discontinuities limit the global
experimental orders, so these are not formal smooth-flow accuracy estimates.

## Strong pressure jump

For the strong-shock trial I used states $(1,0,1000)$ and $(1,0,0.01)$ with
$\gamma=1.4$. I reduced the CFL number to 0.7 and stopped at $t=0.01$, keeping
the waves away from the boundaries. The pressure ratio is $10^5$.

| Diagnostic | Measured value |
|---|---:|
| Density $L_1$ error | $1.089200\times10^{-1}$ |
| Pressure $L_1$ error | $6.699012$ |
| Minimum density | $0.567814$ |
| Minimum pressure | $0.01$ |
| Relative mass boundary-budget residual | $1.414\times10^{-12}$ |
| Relative energy boundary-budget residual | $3.959\times10^{-12}$ |

The method stays positive and captures the wave pattern, but the narrow dense
post-shock shell is broadened and its peak underestimated. The pointwise errors
make this a useful target for reconstruction improvements.

## Translating contact

For the contact trial I used states $(1,1,1)$ and $(0.125,1,1)$, initially
separated at $x=0.3$, and evolved them to $t=0.2$. The exact contact is then at
$x=0.5$.

| Diagnostic | Measured value |
|---|---:|
| Density $L_1$ error | $1.936408\times10^{-2}$ |
| 1--99% density-transition width | $0.135$ |
| Velocity $L_\infty$ error | $1.89\times10^{-15}$ |
| Pressure $L_\infty$ error | $1.11\times10^{-15}$ |
| Relative mass boundary-budget residual | $4.66\times10^{-15}$ |

HLL preserves constant velocity and pressure while strongly diffusing density.
Raw totals change because unequal-density flow crosses the open boundaries;
the boundary-budget residual is the relevant conservation measure.

## Strong rarefaction

For the rarefaction trial I used symmetric states $(1,-2,0.4)$ and
$(1,2,0.4)$ with $\gamma=1.4$, evolved to $t=0.15$ at CFL 0.7. I specifically
checked the central cells because that is where density and pressure approach
their smallest values.

| Diagnostic | Measured value |
|---|---:|
| Density $L_1$ error | $1.060086\times10^{-2}$ |
| Pressure $L_1$ error | $6.266514\times10^{-3}$ |
| Minimum density | $0.016492$ |
| Minimum pressure | $0.004290$ |
| Relative mass boundary-budget residual | $8.14\times10^{-11}$ |
| Relative energy boundary-budget residual | $1.33\times10^{-10}$ |

No density or pressure floor was used. The numerical centre remains positive,
though the first-order method overestimates its exact thermodynamic state and
produces substantial specific-internal-energy diffusion.

## Current validation boundary

No smooth-wave convergence test, second-order reconstruction, alternative
Riemann flux, or multidimensional benchmark has been run. The current evidence
supports the stated first-order HLL baseline only. The next validation target
is a smooth periodic wave alongside MUSCL and RK2.
