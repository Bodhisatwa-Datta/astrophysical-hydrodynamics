# Numerical methods

## Finite-volume discretisation

Cell averages are updated with the conservative formula

$$
\mathbf{U}_i^{n+1}=\mathbf{U}_i^n
-\frac{\Delta t}{\Delta x}
\left(\widehat{\mathbf{F}}_{i+1/2}-
\widehat{\mathbf{F}}_{i-1/2}\right).
$$

The current reconstruction is piecewise constant, so the left and right
states at each face are the adjacent cell averages. Time integration is one
forward-Euler step. This path remains available as the first-order baseline.

## MUSCL reconstruction

The second-order path reconstructs primitive variables
$\mathbf{W}=(\rho,u,p)^T$. For each cell,

$$
\Delta\mathbf{W}_i=
\phi(\mathbf{W}_i-\mathbf{W}_{i-1},
     \mathbf{W}_{i+1}-\mathbf{W}_i),
$$

and the interface states are

$$
\mathbf{W}_{i+1/2}^{L}=\mathbf{W}_i+\frac{1}{2}\Delta\mathbf{W}_i,
\qquad
\mathbf{W}_{i+1/2}^{R}=\mathbf{W}_{i+1}-\frac{1}{2}\Delta\mathbf{W}_{i+1}.
$$

Reconstructing primitive variables makes density and pressure limiting direct
and avoids mixing the thermodynamic constraint into conserved-component
slopes. Reconstructed states are checked for positive density and pressure;
they are not clipped.

The available TVD limiters are minmod, monotonized central (MC), and van Leer.
Minmod is the most dissipative of the three. MC permits a steeper monotone
slope, while van Leer uses a smooth harmonic mean when neighboring slopes have
the same sign. Every limiter returns zero when the one-sided slopes disagree.

## Time integration

MUSCL is paired with the two-stage strong-stability-preserving Runge--Kutta
method

$$
\mathbf{U}^{(1)}=\mathbf{U}^{n}+\Delta t\,L(\mathbf{U}^{n}),
$$

$$
\mathbf{U}^{n+1}=\frac{1}{2}\mathbf{U}^{n}
+\frac{1}{2}\left[\mathbf{U}^{(1)}+\Delta t\,L(\mathbf{U}^{(1)})\right].
$$

Ghost cells and interface fluxes are recomputed at both stages. Using RK2 with
MUSCL keeps temporal and spatial accuracy consistent for smooth flow.

## Approximate Riemann fluxes

### HLL

For left/right states $\mathbf{U}_L$ and $\mathbf{U}_R$, signal speeds are
estimated as

$$
S_L=\min(u_L-c_L,u_R-c_R),\qquad
S_R=\max(u_L+c_L,u_R+c_R).
$$

The HLL flux is

$$
\widehat{\mathbf{F}}_{\rm HLL}=
\begin{cases}
\mathbf{F}_L,&S_L\geq0,\\
\dfrac{S_R\mathbf{F}_L-S_L\mathbf{F}_R+S_LS_R(\mathbf{U}_R-\mathbf{U}_L)}{S_R-S_L},
&S_L<0<S_R,\\
\mathbf{F}_R,&S_R\leq0.
\end{cases}
$$

HLL is robust and captures shocks, but its two-wave model does not explicitly
resolve the contact wave and therefore smears contacts.

### Rusanov

The local Lax--Friedrichs or Rusanov flux is

$$
\widehat{\mathbf{F}}_{\rm Rus}=\frac{1}{2}(\mathbf{F}_L+\mathbf{F}_R)
-\frac{1}{2}a_{\max}(\mathbf{U}_R-\mathbf{U}_L),
$$

where

$$a_{\max}=\max(|u_L|+c_L,|u_R|+c_R).$$

Its single maximum-speed dissipation makes it compact and robust, but usually
more diffusive than HLL or HLLC around contacts and narrow post-shock features.

### HLLC

HLLC restores the missing middle wave. Using the same Davis outer speeds as
HLL, the contact speed is

$$
S_M=\frac{p_R-p_L+\rho_Lu_L(S_L-u_L)-\rho_Ru_R(S_R-u_R)}
{\rho_L(S_L-u_L)-\rho_R(S_R-u_R)}.
$$

Left and right star states are constructed across $S_L$, $S_M$, and $S_R$,
and their fluxes follow from the Rankine--Hugoniot relation

$$\mathbf{F}_{K}^{*}=\mathbf{F}_{K}+S_K(\mathbf{U}_{K}^{*}-\mathbf{U}_{K}).$$

The implementation computes star pressure from both sides and averages the two
algebraically identical values to reduce floating-point asymmetry. HLLC exactly
recovers the physical flux of a stationary constant-pressure contact. It costs
more arithmetic and is not assumed to be more robust or accurate for every
wave pattern.

## Stability and boundaries

The timestep is recomputed before every update:

$$
\Delta t=C_{\rm CFL}\frac{\Delta x}{\max_i(|u_i|+c_{s,i})}.
$$

Only the final step is shortened to land exactly on the requested output time.
The default CFL number is 0.8. Transmissive boundaries copy the nearest active
cell into each ghost cell, corresponding to a zero-gradient extrapolation.
Periodic boundaries copy active cells from the opposite edge and are used for
the smooth entropy-wave experiment.

For the current open-boundary Riemann problems, conservation is assessed using
the finite-volume boundary budget

$$
\Delta\sum_i\mathbf{U}_i\Delta x
=\left(\mathbf{F}_{\rm left}-\mathbf{F}_{\rm right}\right)\Delta t.
$$

The reported residual subtracts this expected flux contribution. This simple
far-field expression is used only while the initial left and right states still
occupy the domain boundaries.

## Exact validation solution

The internal exact Riemann solver iterates for the star-region pressure using
the shock and rarefaction pressure functions, then samples the self-similar
solution in $\xi=(x-x_0)/t$. Its construction follows the standard method in
E. F. Toro, *Riemann Solvers and Numerical Methods for Fluid Dynamics*. It is
used only for comparisons and never to evolve the numerical state.

## Error norms and observed order

For a cell-centred numerical field $q_i$ and exact value $q_i^{\rm exact}$,
the reported discrete integral error is

$$E_1=\Delta x\sum_i|q_i-q_i^{\rm exact}|.$$

Between two grid spacings, the experimental order is

$$p=\frac{\log(E_h/E_{h/2})}{\log 2}.$$

Shock-tube solutions contain discontinuities, so their global observed order
need not equal the method's formal order on smooth solutions. Smooth periodic
entropy waves are therefore measured separately in one and two dimensions.

## Two-dimensional update

The current 2D solver uses an unsplit conservative update

$$
\mathbf{U}_{i,j}^{n+1}=\mathbf{U}_{i,j}^{n}
-\frac{\Delta t}{\Delta x}
(\widehat{\mathbf{F}}_{i+1/2,j}-\widehat{\mathbf{F}}_{i-1/2,j})
-\frac{\Delta t}{\Delta y}
(\widehat{\mathbf{G}}_{i,j+1/2}-\widehat{\mathbf{G}}_{i,j-1/2}).
$$

The baseline uses piecewise-constant interface states and forward Euler. The
second-order path applies the same primitive-variable TVD reconstruction
independently along x and y, then combines both flux divergences in a single
method-of-lines operator. SSP-RK2 evaluates that complete unsplit operator at
both stages. This is direction-by-direction reconstruction, not dimensional
time splitting; it does not include a corner-transport or transverse predictor.

For y faces, normal and tangential momentum are rotated into the same ordering
used at x faces, the chosen HLL, HLLC, or Rusanov construction is applied, and
the momentum fluxes are rotated back. HLLC carries tangential velocity into its
star states.

The multidimensional timestep is

$$
\Delta t=C_{\rm CFL}\left[
\max_{i,j}\left(
\frac{|v_x|+c_s}{\Delta x}+\frac{|v_y|+c_s}{\Delta y}
\right)\right]^{-1}.
$$

Outflow boundaries use zero-gradient extrapolation. Periodic boundaries copy
the opposite active cells. Reflective boundaries reverse only the momentum
normal to the wall while copying density, tangential momentum, and energy.

## Kelvin--Helmholtz benchmark and diagnostics

The benchmark uses a periodic unit square with two smooth shear layers. A
window function

$$
L(y)=\frac{1}{2}\left[
\tanh\left(\frac{y-0.25}{a}\right)-
\tanh\left(\frac{y-0.75}{a}\right)
\right],\qquad a=0.025,
$$

sets

$$
\rho=1+L,\qquad v_x=-0.5+L,\qquad p=2.5.
$$

The transverse seed is

$$
v_y=0.01\sin(4\pi x)\left[
e^{-d_{0.25}(y)^2/(2\sigma^2)}+
e^{-d_{0.75}(y)^2/(2\sigma^2)}
\right],\qquad \sigma=0.05,
$$

where each $d_{y_0}$ is the shortest signed periodic distance from $y_0$.
This selects two wavelengths across the domain and avoids a grid-scale jump in
density or horizontal velocity.

Cell-centred vorticity is calculated with periodic centred differences,

$$
\omega_z=\frac{\partial v_y}{\partial x}-
\frac{\partial v_x}{\partial y}.
$$

The primary growth diagnostic is transverse kinetic energy,

$$
K_y=\sum_{i,j}\frac{1}{2}\rho_{i,j}v_{y,i,j}^2\,\Delta x\Delta y.
$$

RMS vertical velocity, RMS and maximum vorticity, density standard deviation,
thermodynamic minima, conserved totals, runtime, and step count are recorded at
fixed output times. Since no physical viscosity is present, resolution changes
alter the effective numerical dissipation and must be considered when
interpreting growth and small-scale structure.

## Constant gravity and hydrostatic walls

Uniform gravity is added to the same method-of-lines operator as the flux
divergence, so both SSP-RK2 stages recompute the momentum and energy sources.
Besides the acoustic CFL limit, the timestep is restricted by

$$
\Delta t_g=\sqrt{C_{\rm CFL}\frac{\min(\Delta x,\Delta y)}{|\mathbf g|}}.
$$

At a hydrostatic horizontal wall, density and tangential momentum are reflected
and normal momentum changes sign. Ghost-cell pressure is extrapolated using

$$
p_{j+1}-p_j=\frac{1}{2}(\rho_j+\rho_{j+1})g_y\Delta y,
$$

which is consistent with $dp/dy=\rho g_y$. This reduces boundary-driven motion
but does not make the complete finite-volume method exactly well balanced;
unperturbed control calculations are therefore retained and measured.

## Rayleigh--Taylor benchmark

The unit-square benchmark uses $g_y=-0.5$, a light density of 1 below a heavy
density of 2, and a smooth transition of width $a=0.025$ at $y_0=0.5$:

$$
\rho(y)=\bar\rho+\Delta\rho\tanh\left(\frac{y-y_0}{a}\right),
$$

where $\bar\rho=1.5$ and $\Delta\rho=0.5$. Integrating hydrostatic balance
analytically gives

$$
p(y)=p_0+g_y\left[
\bar\rho(y-y_0)+\Delta\rho\,a
\ln\cosh\left(\frac{y-y_0}{a}\right)
\right],\qquad p_0=2.5.
$$

The perturbed calculation starts with

$$
v_y=0.0025\sin(4\pi x)
\exp\left[-\frac{(y-y_0)^2}{2(0.05)^2}\right],
$$

while the control uses exactly zero velocity. The density-$1.5$ contour is
interpolated in every column; its maximum and minimum provide bubble height and
spike depth. An exponential rate is fitted to the mean interface amplitude
between $t=0.8$ and $t=2.2$. Gas-plus-potential energy, vertical kinetic energy,
vorticity, conservation, and the unperturbed RMS vertical velocity are also
recorded.

## Sedov--Taylor energy deposition and shock tracking

The blast uses a uniform ambient state $(\rho,p)=(1,10^{-5})$ on
$[-0.5,0.5]^2$. One unit of thermal energy is deposited inside radius
$r_{\rm inj}=0.05$ with the compact kernel

$$
w(r)=
\begin{cases}
[1-(r/r_{\rm inj})^2]^2,&r<r_{\rm inj},\\
0,&r\ge r_{\rm inj}.
\end{cases}
$$

On each grid, the sampled weights are normalized by
$\sum_{i,j}w_{i,j}\Delta x\Delta y$. The discrete injected energy is therefore
exactly $E_0=1$ rather than varying with the number of cells inside the
injection region. The energy is thermal, velocities initially vanish, and no
single cell receives the entire explosion.

Cell values are averaged in annuli of width $\min(\Delta x,\Delta y)$. The
shock radius is the location of the steepest outward decrease in the radial
density profile, with a local quadratic interpolation of the gradient minimum.
A power law $R=Ct^\alpha$ is fitted from $t=0.015$ to $0.05$, after the shock
has expanded well beyond the injection radius. The measured $\alpha$ is
compared with the two-dimensional similarity prediction $1/2$.

Angular symmetry is measured independently by repeating the radial-gradient
location in 16 equal polar sectors. The standard deviation and peak-to-peak
range of those sector radii are normalized by their mean. These measurements
include hydrodynamic anisotropy and finite sampling error from the detector.
