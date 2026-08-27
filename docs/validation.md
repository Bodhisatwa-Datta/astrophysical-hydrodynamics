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
rarefaction corner, contact discontinuity, and shock, and this baseline
configuration is first order.

## Second-order 2D extension

I extended primitive-variable MUSCL reconstruction separately along both
coordinate directions and paired the complete unsplit flux-divergence operator
with SSP-RK2. I first repeated the uniform periodic calculation using
MUSCL-MC/RK2. After $t=0.2$, the maximum change in any conserved cell value and
the reported changes in all four integrated conserved quantities were zero.

For a smooth genuinely two-dimensional measurement, I used

$$
\rho(x,y,0)=1+0.2\sin[2\pi(x+y)],\qquad
(v_x,v_y)=(0.7,0.3),\qquad p=1.
$$

This is an exact entropy wave. On the periodic unit square it returns to its
initial profile at $t=1$. I ran the first-order configuration and all three
MUSCL limiters at $16^2$, $32^2$, and $64^2$, always with HLL and CFL 0.4.

| Method | $L_1$ at $16^2$ | $L_1$ at $64^2$ | Order, $32^2$ to $64^2$ |
|---|---:|---:|---:|
| First order + Euler | $1.180130\times10^{-1}$ | $6.401181\times10^{-2}$ | 0.578 |
| MUSCL minmod + RK2 | $6.212038\times10^{-2}$ | $8.852759\times10^{-3}$ | 1.146 |
| MUSCL MC + RK2 | $1.759974\times10^{-2}$ | $1.202684\times10^{-3}$ | 1.698 |
| MUSCL van Leer + RK2 | $2.992421\times10^{-2}$ | $2.566743\times10^{-3}$ | 1.814 |

MC and van Leer converge much faster than the first-order baseline and trend
toward second order. I do not treat three relatively modest resolutions as a
perfect asymptotic proof, especially because TVD limiting near smooth extrema
reduces the measured rate. Across these runs, the largest relative change in
mass, either momentum component, or energy was $3.7\times10^{-16}$. Pressure
remained unity to roundoff and density remained positive without a floor.

I also repeated the rotated Sod calculation with MUSCL-MC/RK2:

| Grid | Density $L_1$, x-normal | Density $L_1$, y-normal | Order | Rotational $L_\infty$ difference |
|---:|---:|---:|---:|---:|
| $32^2$ | $1.222166\times10^{-2}$ | $1.222166\times10^{-2}$ | -- | $0.0$ |
| $64^2$ | $7.459356\times10^{-3}$ | $7.459356\times10^{-3}$ | 0.712 | $0.0$ |
| $128^2$ | $4.075024\times10^{-3}$ | $4.075024\times10^{-3}$ | 0.872 | $0.0$ |

The shock-tube order remains below the smooth-wave rate because the exact
solution contains discontinuities. At every resolution, density and pressure
stayed at or above 0.125 and 0.1, and the rotated solutions agreed exactly
after transposition and velocity-component exchange.

## Smooth Kelvin--Helmholtz shear layer

I initialized two periodic hyperbolic-tangent shear layers of width 0.025. The
central fluid has density 2 and horizontal velocity $+0.5$; the outer fluid has
density 1 and velocity $-0.5$. Pressure is uniformly 2.5. I seeded both layers
with a Gaussian-localised transverse perturbation
$0.01\sin(4\pi x)$, then evolved the flow to $t=1.5$ with HLLC,
MUSCL-MC/RK2, CFL 0.35, and periodic boundaries in both directions.

I repeated the calculation at $64^2$, $96^2$, and $128^2$. The initial
transverse kinetic energy was $6.646702\times10^{-6}$ at every resolution.

| Grid | Steps | $K_y(t=1.5)$ | Growth factor | RMS $v_y$ | Maximum $|\omega_z|$ | Runtime |
|---:|---:|---:|---:|---:|---:|---:|
| $64^2$ | 1188 | $3.973554\times10^{-3}$ | 597.8 | 0.07583 | 18.246 | 14.8 s |
| $96^2$ | 1787 | $5.839117\times10^{-3}$ | 878.5 | 0.09254 | 19.881 | 31.3 s |
| $128^2$ | 2373 | $6.564919\times10^{-3}$ | 987.7 | 0.09838 | 21.915 | 81.5 s |

The transverse energy first falls during adjustment of the imposed
perturbation, then grows by nearly three orders of magnitude. By $t=1.5$, two
paired billows deform each shear layer, as selected by the seeded wavelength.
The $96^2$ and $128^2$ RMS vertical velocities differ by about 6%, and their
transverse energies differ by about 11%. The $64^2$ calculation is more
diffusive and has 39% less transverse energy than the $128^2$ calculation.
This is clear resolution dependence rather than evidence of full convergence.

The largest absolute mass drift was $2.22\times10^{-16}$ and the largest
absolute energy drift was $8.88\times10^{-16}$. Total horizontal momentum was
unchanged to reported precision; total vertical momentum, whose exact initial
value is zero, drifted by at most $5.10\times10^{-18}$. At $128^2$, the minimum
density and pressure were 0.96296 and 2.16754, so the run remained positive
without floors. The maximum vorticity continues to rise with resolution,
which is expected when an inviscid calculation resolves a thinner shear layer.

I have not compared the growth curve with an independently derived linear
eigenmode for this finite-width, compressible, double-layer setup. I therefore
interpret the calculation as a reproducible nonlinear instability and
resolution study, not a measurement of an analytic growth rate. With no
physical viscosity, the dissipation visible here is numerical.

The second-order path has now been measured on periodic entropy waves in one
and two dimensions, four 1D Riemann problems, rotated 2D shock tubes, and a
smooth Kelvin--Helmholtz instability. The instability grows reproducibly and
the two finest runs are approaching one another in global diagnostics, but the
vorticity field and growth amplitude are not grid-converged.

## Rayleigh--Taylor instability and hydrostatic control

I added constant gravity with $g_y=-0.5$ and initialized a smooth transition
from density 1 below to density 2 above. Pressure was obtained by analytically
integrating $dp/dy=\rho g_y$, with $p=2.5$ at $y=0.5$. Horizontal boundaries
are periodic and the vertical walls use reflective velocity with hydrostatic
pressure extrapolation. The perturbed run begins with a localized
$0.0025\sin(4\pi x)$ vertical velocity; its paired control begins with exactly
zero velocity.

I evolved both cases to $t=2.5$ using HLLC, MUSCL-MC/RK2, and CFL 0.3 at three
resolutions. The perturbed results were:

| Grid | Bubble height | Spike depth | $K_y$ | RMS $v_y$ | Fitted interface rate |
|---:|---:|---:|---:|---:|---:|
| $64^2$ | 0.01244 | 0.01268 | $1.1592\times10^{-5}$ | 0.003978 | 1.410 |
| $96^2$ | 0.01264 | 0.01289 | $1.1994\times10^{-5}$ | 0.004044 | 1.409 |
| $128^2$ | 0.01282 | 0.01310 | $1.2088\times10^{-5}$ | 0.004059 | 1.397 |

I fitted the mean of bubble height and spike depth to an exponential over
$0.8\le t\le2.2$. The three measured rates agree within 0.9%. For context, the
sharp-interface incompressible estimate
$\sqrt{A|g|k}=1.447$, using $A=1/3$ and $k=4\pi$, is 2.6--3.5% larger. This is
only a reference scale: the simulated interface has finite width, the gas is
compressible, and the imposed velocity is not an exact eigenmode.

The two finest bubble heights differ by 1.4%, their spike depths by 1.6%, and
their vertical kinetic energies by 0.8%. This is much closer agreement than in
the Kelvin--Helmholtz study, although the calculation is still in the linear
regime and does not establish convergence of nonlinear mixing.

The unperturbed controls quantify residual hydrostatic imbalance:

| Grid | RMS $v_y$ at $t=2.5$ | Maximum $|v_y|$ | Relative gas+potential energy drift |
|---:|---:|---:|---:|
| $64^2$ | $1.395\times10^{-5}$ | $5.688\times10^{-5}$ | $1.277\times10^{-7}$ |
| $96^2$ | $4.413\times10^{-6}$ | $1.656\times10^{-5}$ | $3.913\times10^{-8}$ |
| $128^2$ | $2.359\times10^{-6}$ | $7.650\times10^{-6}$ | $1.670\times10^{-8}$ |

The control drift decreases strongly with refinement and is over three orders
of magnitude smaller than the perturbed RMS velocity on the finest grid. The
perturbed gas-plus-potential-energy drift remained below $5.5\times10^{-9}$ in
magnitude. All runs retained positive states; the lowest recorded density and
pressure were 0.99921 and 2.00719.

## Two-dimensional Sedov--Taylor blast

I initialized a uniform ambient medium with density 1, pressure $10^{-5}$, and
$\gamma=1.4$ on $[-0.5,0.5]^2$. I deposited one unit of thermal energy through
a compact polynomial kernel inside radius 0.05. I normalized the sampled
kernel separately on each grid, so the excess energy above the ambient value
was exactly 1 at $64^2$, $96^2$, and $128^2$. The runs used MUSCL-MC/RK2, HLL,
CFL 0.25, and outflow boundaries; the shock remained well inside the domain.

I measured the shock from the steepest outward density decrease in annular
profiles and fitted $R=Ct^\alpha$ over $0.015\le t\le0.05$. I independently
evaluated the closed-form cylindrical similarity solution for the same
$E_0$, $\rho_0$, and $\gamma$. Its dimensionless coefficient is
$\xi_2=1.00403$, giving $R_s(0.05)=0.22451$.

| Grid | Steps | $R(0.05)$ | Radius error | Fitted $\alpha$ | Error from $1/2$ | Angular scatter |
|---:|---:|---:|---:|---:|---:|---:|
| $64^2$ | 239 | 0.23542 | $+4.86\%$ | 0.4630 | $-7.40\%$ | $1.61\%$ |
| $96^2$ | 375 | 0.23180 | $+3.25\%$ | 0.4727 | $-5.45\%$ | $0.97\%$ |
| $128^2$ | 510 | 0.22958 | $+2.26\%$ | 0.4826 | $-3.49\%$ | $0.80\%$ |

The exponent moves monotonically toward the two-dimensional similarity value
$1/2$. I did not use the familiar spherical $2/5$ exponent because it belongs
to a three-dimensional blast. The two finest final radii differ by 0.96%.
Their radial profiles place the density shell and pressure drop at nearly the
same radius, while the increasing density peak shows that the thin swept-up
shell is still resolution sensitive.

I compared the annular means with the exact density, radial velocity, and
pressure profiles over $0\le r\le1.25R_s$. The relative, radius-weighted $L_1$
errors decrease consistently:

| Grid | Density | Radial velocity | Pressure |
|---:|---:|---:|---:|
| $64^2$ | 23.38% | 21.43% | 18.14% |
| $96^2$ | 20.79% | 9.14% | 12.26% |
| $128^2$ | 14.90% | 8.44% | 8.89% |

The density error converges more slowly because the exact solution jumps from
six times the ambient density immediately behind the shock to the undisturbed
state, whereas the numerical method spreads that shell over several cells.

For the angular measurement I located the shock independently in 16 sectors.
Its relative standard deviation decreases by a factor of two from $64^2$ to
$128^2$; the finest peak-to-peak angular range is 1.66% of the mean radius.
The residual fourfold imprint visible on the coarse Cartesian grid weakens
under refinement but has not vanished.

Mass changes were zero to reported precision. The largest relative total-
energy change was $2.23\times10^{-16}$, and absolute momentum drift stayed
below $1.74\times10^{-17}$. The minimum density reached 0.03760 in the central
rarefaction and the minimum pressure remained at the positive ambient value
$10^{-5}$; no floor or clipping was applied.

I used the standard strong-shock, zero-ambient-pressure similarity solution for
the reference profiles while retaining $p_0=10^{-5}$ outside the analytic
shock for plotting. The numerical run begins from a finite injection region,
not an ideal point explosion, so early-time and finite-resolution offsets are
expected. The closed-form equations and normalization follow J. R. Kamm,
[Evaluation of the Sedov-von Neumann-Taylor Blast Wave
Solution](https://cococubed.com/papers/kamm_2000.pdf), LA-UR-00-6055 (2000).

## Current validation boundary

The planned sequence now includes validated 1D shock tubes, reconstruction and
flux comparisons, smooth-order studies, a second-order 2D core,
Kelvin--Helmholtz growth, controlled Rayleigh--Taylor growth, and a
two-dimensional Sedov blast. Remaining limitations are quantitative rather
than hidden: Kelvin--Helmholtz is not grid-converged, gravity is not exactly
well balanced, Rayleigh--Taylor has not reached nonlinear mixing, and Sedov has
a broadened density shell and a residual 2.26% shock-radius offset on the finest
grid.
