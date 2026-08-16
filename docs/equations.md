# Governing equations

## One-dimensional Euler system

The implemented system is

$$
\partial_t \mathbf{U} + \partial_x \mathbf{F}(\mathbf{U})=0,
$$

with conserved state and physical flux

$$
\mathbf{U}=
\begin{pmatrix}\rho\\ \rho u\\ E\end{pmatrix},
\qquad
\mathbf{F}=
\begin{pmatrix}\rho u\\ \rho u^2+p\\ u(E+p)\end{pmatrix}.
$$

Here $\rho>0$ is mass density, $u$ is velocity, $p>0$ is thermal
pressure, and $E$ is total energy density. The ideal-gas closure is

$$
E=\frac{p}{\gamma-1}+\frac{1}{2}\rho u^2,
\qquad
p=(\gamma-1)\left(E-\frac{1}{2}\rho u^2\right),
$$

and the adiabatic sound speed is

$$c_s=\sqrt{\gamma p/\rho}.$$

The adiabatic index is configurable and must exceed one. The Sod benchmark
uses $\gamma=1.4$.

## Admissible states

The code requires finite values, strictly positive density, and strictly
positive pressure. It does not silently impose floors. Primitive-to-conserved
and conserved-to-primitive conversion both enforce these conditions so that a
loss of physical admissibility is exposed at its point of detection.

## Conservation diagnostics

On a uniform grid the discrete totals are

$$M=\sum_i\rho_i\Delta x,\qquad
P=\sum_i(\rho u)_i\Delta x,\qquad
E_{\rm tot}=\sum_iE_i\Delta x.$$

For an open domain these totals can change through boundary fluxes. In the Sod
run reported here the waves have not reached the boundaries at $t=0.2$;
momentum nevertheless changes because unequal far-field pressures exert a net
flux across the two domain boundaries. Mass and energy boundary fluxes remain
zero for the stationary far-field states.

