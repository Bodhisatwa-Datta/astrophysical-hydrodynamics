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

## Smooth entropy-wave convergence

To measure formal accuracy away from discontinuities, I advected

$$
\rho(x,0)=1+0.2\sin(2\pi x),\qquad u=1,\qquad p=1
$$

through one periodic unit domain. Constant velocity and pressure make this an
exact Euler solution, and after $t=1$ the density profile should return to its
initial position. I used CFL 0.4 for every method and repeated the calculation
at 32, 64, 128, and 256 cells.

| Method | $L_1$ at 32 cells | $L_1$ at 256 cells | Order, 128 to 256 |
|---|---:|---:|---:|
| First order + Euler | $5.944591\times10^{-2}$ | $9.614786\times10^{-3}$ | 0.944 |
| MUSCL minmod + RK2 | $1.297335\times10^{-2}$ | $3.698817\times10^{-4}$ | 1.864 |
| MUSCL MC + RK2 | $4.030144\times10^{-3}$ | $8.687254\times10^{-5}$ | 1.969 |
| MUSCL van Leer + RK2 | $6.843088\times10^{-3}$ | $1.173868\times10^{-4}$ | 2.038 |

The first-order result approaches order one, while MC and van Leer approach
order two on the finest pair. Minmod also trends toward second order but loses
more accuracy near the smooth extrema. Periodic mass changes remained at or
below $4.5\times10^{-16}$ in the recorded runs.

## Reconstruction comparison on discontinuous flow

I then reran all four Riemann problems at 400 cells and CFL 0.4 so the
first-order and MUSCL results were compared under the same timestep policy.
The table gives density $L_1$ error.

| Method | Sod | Contact | Strong shock | Rarefaction |
|---|---:|---:|---:|---:|
| First order | $7.811995\times10^{-3}$ | $1.990704\times10^{-2}$ | $1.185164\times10^{-1}$ | $1.511274\times10^{-2}$ |
| MUSCL minmod | $2.467526\times10^{-3}$ | $7.728979\times10^{-3}$ | $6.129825\times10^{-2}$ | $4.545778\times10^{-3}$ |
| MUSCL MC | $1.496039\times10^{-3}$ | $4.359090\times10^{-3}$ | $3.981569\times10^{-2}$ | $1.677545\times10^{-3}$ |
| MUSCL van Leer | $1.664298\times10^{-3}$ | $5.025405\times10^{-3}$ | $4.528675\times10^{-2}$ | $2.342627\times10^{-3}$ |

The 1--99% contact width fell from 0.140 for the first-order run to 0.060 with
minmod, 0.0275 with MC, and 0.035 with van Leer. All recorded densities and
pressures remained positive without floors. MC gave the lowest errors in this
particular comparison, but that is not evidence that it is universally best;
the limiter choice still trades smoothness, compression, and robustness.

## HLL, HLLC, and Rusanov comparison

I compared the three fluxes at 400 cells and CFL 0.4 without changing the
initial conditions or output times. Each flux was run once with the
piecewise-constant/forward-Euler configuration and once with MUSCL-MC/RK2.
Complete per-run errors, minima, conservation residuals, step counts, and
runtimes are stored in `benchmarks/riemann_solvers/comparison.csv`.

### Density error with the first-order configuration

| Flux | Sod | Contact | Strong shock | Rarefaction |
|---|---:|---:|---:|---:|
| HLL | $7.811995\times10^{-3}$ | $1.990704\times10^{-2}$ | $1.185164\times10^{-1}$ | $1.511274\times10^{-2}$ |
| HLLC | $7.376189\times10^{-3}$ | $1.486001\times10^{-2}$ | $1.167735\times10^{-1}$ | $1.517544\times10^{-2}$ |
| Rusanov | $1.113357\times10^{-2}$ | $2.524677\times10^{-2}$ | $1.443782\times10^{-1}$ | $1.442088\times10^{-2}$ |

### Density error with MUSCL-MC/RK2

| Flux | Sod | Contact | Strong shock | Rarefaction |
|---|---:|---:|---:|---:|
| HLL | $1.496039\times10^{-3}$ | $4.359090\times10^{-3}$ | $3.981569\times10^{-2}$ | $1.677545\times10^{-3}$ |
| HLLC | $1.446073\times10^{-3}$ | $3.853941\times10^{-3}$ | $3.992733\times10^{-2}$ | $1.677466\times10^{-3}$ |
| Rusanov | $1.872238\times10^{-3}$ | $4.819122\times10^{-3}$ | $4.469082\times10^{-2}$ | $1.744289\times10^{-3}$ |

For the translating contact, the first-order 1--99% widths were 0.140 (HLL),
0.0975 (HLLC), and 0.175 (Rusanov). With MUSCL-MC/RK2 they narrowed to 0.0275,
0.0250, and 0.0275. HLLC therefore gave the clearest contact advantage, while
the reconstruction reduced the difference between fluxes considerably.

All 24 runs retained positive density and pressure without floors. On the
strong shock, HLL and HLLC were effectively tied in density error with the
second-order configuration; HLLC had the lower pressure error, while Rusanov
remained more diffusive. In the first-order rarefaction, Rusanov happened to
give the smallest density error but not the smallest pressure error. I therefore
do not treat the ranking from one scalar metric as universal.

The mean wall times over the four problems were:

| Configuration | HLL | HLLC | Rusanov |
|---|---:|---:|---:|
| First order | 0.335 s | 0.366 s | 0.345 s |
| MUSCL-MC/RK2 | 0.862 s | 0.964 s | 0.910 s |

In this single local run, HLLC cost roughly 9--12% more than HLL. These timings
are descriptive rather than a hardware-independent performance result.

## First-order 2D Euler core

I began the 2D validation with a uniform periodic state on a $48\times36$ grid:

$$
(\rho,v_x,v_y,p)=(1.1,0.3,-0.2,0.9).
$$

After evolving to $t=0.2$ with HLLC, the largest change in any conserved cell
value was reported as zero, as were the relative changes in total mass, both
momentum components, and total energy. I separately checked the directional
flux identities for HLL, HLLC, and Rusanov and inspected reflective ghost cells
to confirm that only wall-normal momentum changes sign.

I then embedded the Sod problem in a square 2D grid, first with the initial
discontinuity normal to x and then with the entire problem rotated by 90
degrees. Both calculations used first-order unsplit HLL, CFL 0.4, and final
time $t=0.1$.

| Grid | Density $L_1$, x-normal | Density $L_1$, y-normal | Order | Rotational $L_\infty$ difference |
|---:|---:|---:|---:|---:|
| $32^2$ | $2.418633\times10^{-2}$ | $2.418633\times10^{-2}$ | -- | $0.0$ |
| $64^2$ | $1.875955\times10^{-2}$ | $1.875955\times10^{-2}$ | 0.367 | $0.0$ |
| $128^2$ | $1.322181\times10^{-2}$ | $1.322181\times10^{-2}$ | 0.505 | $0.0$ |

At $64^2$ and $128^2$, mass and total-energy changes were zero to reported
precision. The small $32^2$ changes, $2.56\times10^{-11}$ in mass and
$4.22\times10^{-11}$ in energy, are consistent with the long numerical tails
of the coarse first-order solution reaching the outflow boundary. Density and
pressure remained positive without floors at every recorded resolution.

The identical rotated profiles show that the directional momentum mapping and
flux indexing are symmetric for this grid-aligned problem. The low observed
orders are not a smooth-flow accuracy measurement: the solution contains a
rarefaction corner, contact discontinuity, and shock, and the 2D method is
currently first order.

## Current validation boundary

The second-order result has been demonstrated for one smooth entropy wave and
the limiter and flux behavior has been measured on four 1D Riemann problems.
The 2D homogeneous Euler core has passed uniform-flow and grid-aligned rotated
shock-tube checks, but it remains first order. The next target is an unsplit 2D
MUSCL-RK2 extension and smooth-advection convergence study before any
Kelvin--Helmholtz result is attempted.
