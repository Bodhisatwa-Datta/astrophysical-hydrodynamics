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
forward-Euler step. The resulting method is first order in smooth regions.

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

## Exact validation solution

The internal exact Riemann solver iterates for the star-region pressure using
the shock and rarefaction pressure functions, then samples the self-similar
solution in $\xi=(x-x_0)/t$. Its construction follows the standard method in
E. F. Toro, *Riemann Solvers and Numerical Methods for Fluid Dynamics*. It is
used only for comparisons and never to evolve the numerical state.

