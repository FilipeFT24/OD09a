To illustrate the **dual consistent** discretisation of the source term, we consider the following nonlinear *elliptic* ODE proposed by Oliver *et al.* (see "About" section for further reference):

$$
-((\nu+u)u_{x})_{x}-cu_{x}u_{x}=g(x)\quad\text{for }x\in[0,1]\text{ and } u(0)=u(1)=0,
$$

where $\nu=1$ and $c=1/2$.

The source (or forcing term) $g(x)$ is chosen to be

$$
g(x)=\pi^{2}[(\nu+\sin(\pi x))\sin(\pi x)-(1+c)\cos^{2}(\pi x)]
$$

to obtain the exact solution

$$
u(x)=\sin(\pi x).
$$

<br />

This case is reproduced here (in this repository) using the BR2 method **with** and **without** lifting contribution.

**NOTE 1):** We use the Newton (or Newton-Raphson) method to solve $J(u^{k})\delta_{u}=R(u^{k})$, with $\delta_{u}=u^{k}-u^{k+1}$ until $\lVert\delta_{u}\rVert_{\infty}<10^{-12}$ (or some other user-defined tolerance).

**NOTE 2):** The construction of the spatial operators follows a similar strategy to that employed in the uhainafil-toy1d repository: https://github.com/FilipeFT24/uhainafil-toy1d **(currently private)**.

**NOTE 3):** There seems to be a bug for the adjoint problem with non-homogeneous boundary conditions! (problem for another day)