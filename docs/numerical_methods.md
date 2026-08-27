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

## HLL interface flux

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
need not equal the method's formal order on smooth solutions. A smooth-wave
test will be required alongside the future second-order scheme.
