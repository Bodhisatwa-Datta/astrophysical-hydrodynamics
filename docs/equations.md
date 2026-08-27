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

## Two-dimensional Euler system

The 2D extension solves

$$
\partial_t\mathbf{U}+\partial_x\mathbf{F}+\partial_y\mathbf{G}=0,
$$

with

$$
\mathbf{U}=
\begin{pmatrix}\rho\\ \rho v_x\\ \rho v_y\\ E\end{pmatrix},
$$

$$
\mathbf{F}=
\begin{pmatrix}
\rho v_x\\ \rho v_x^2+p\\ \rho v_xv_y\\ v_x(E+p)
\end{pmatrix},
\qquad
\mathbf{G}=
\begin{pmatrix}
\rho v_y\\ \rho v_xv_y\\ \rho v_y^2+p\\ v_y(E+p)
\end{pmatrix}.
$$

The total energy is

$$
E=\frac{p}{\gamma-1}+\frac{1}{2}\rho(v_x^2+v_y^2).
$$

With a prescribed constant gravitational acceleration
$\mathbf{g}=(g_x,g_y)$, the implemented source is

$$
\mathbf{S}_{\rm grav}=
\begin{pmatrix}
0\\ \rho g_x\\ \rho g_y\\
\rho v_xg_x+\rho v_yg_y
\end{pmatrix}.
$$

The last component is gravitational work on the gas. For the time-independent
potential $\Phi=-(g_xx+g_yy)$, a closed system should conserve the sum of gas
and potential energy,

$$E_{\rm gas+grav}=\int(E+\rho\Phi)\,dV.$$

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

## Sedov--Taylor dimensional scaling

For an impulsive energy $E_0$ released into uniform density $\rho_0$, dimensional
analysis in $d$ spatial dimensions gives

$$
R(t)=\xi_d(\gamma)
\left(\frac{E_0t^2}{\rho_0}\right)^{1/(d+2)}.
$$

The implemented calculation is two-dimensional Cartesian flow with circular
symmetry, so its expected time exponent is

$$R(t)\propto t^{1/2}.$$

The three-dimensional spherical exponent $2/5$ is not applicable to this
geometry. The coefficient $\xi_2(\gamma)$ depends on the similarity solution
and is not assumed in the current comparison.
