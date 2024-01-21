# Bio Condensates

Energy Functional:

$$
\begin{equation}
    F = \int \mathrm{d}x \left[
        \frac{\varepsilon}{2} |\nabla \phi|^2 +
        \frac{1}{\varepsilon} W(\phi)
    \right]
\end{equation}
$$

with

$$
\begin{equation}
    W(\phi) = 18 \phi^2 \left( 1 - \phi \right)^2
\end{equation}
$$

The double well function $W(\phi)$. Defines what values the $\phi$ can take,
as it has minima at $0, 1$. (Other choices of double well function are possible.)

$L_2$ gradient flow converges to mean curvature flow for
$\varepsilon \rightarrow 0$:

$$
\begin{gathered}
    \varepsilon^{???} \partial_t \phi = - \frac{\delta F}{\delta \phi} \\
    \varepsilon \rightarrow 0: \qquad
    V = -H
\end{gathered}
$$

With $V$ velocity and $H$ mean curvature. Example for Allen Cahn $L_2$ flow of
shrinking circle: Should converge to $\partial_t R = - \frac{1}{R}$ with
R being the circle radius. <br>
This also means that we can choose $\frac{\delta F}{\delta \phi}$ as an approximation
of the curvature.

You can replace $(\nabla \phi)^2$ in the flow equations with
$\frac{1}{\varepsilon} W(\phi)$ since the idea is that both just give
contributions on the interface. However with the double well function we
have numerically more options and scaling with the $\varepsilon$.
With $(\nabla \phi)^2$ you can't really create a weak version, but have
to resort to a simple $\langle (\nabla \phi)^2, \Phi \rangle$ with $\Phi$
being the test function. With $W(\phi)$ you can do stuff like a Taylor
expansion around $x_0 = \phi_n$ and then plug $\phi_{n+1}$ in when plugging into
flow.

# CH Curvature

From [Link to seminar video](https://www.newton.ac.uk/seminar/40850/)

$$
\begin{equation}
    \partial_t \phi = M \nabla^2 \frac{\delta F}{\delta \phi} + R(\phi) + k (\nabla \phi)^4 H(\phi)
\end{equation}
$$

Where $H(\phi)$ is related to the curvature. They chose:

$$
\begin{equation}
    - \nabla \cdot \frac{\nabla \phi}{|\nabla \phi|}
\end{equation}
$$

We choose $H(\phi) = \mu = \frac{\delta F}{\delta \phi}$ for reason detailed above. <br>
$(\nabla \phi)^4$ will be replaced with $W(\phi)$, since its only job is to give contributions on the interface, and $W(\phi)$ does exactly that but is numerically friendly. <br>
We also ignore the $R(\phi)$ part, as they use it to model reactions between two species of droplets. This results in:

$$
\begin{equation}
    \partial_t \phi = M \nabla^2 \frac{\delta F}{\delta \phi} + \frac{k}{\varepsilon} W(\phi) \frac{\delta F}{\delta \phi}
\end{equation}
$$

# Model b+

Original equation:

$$
\begin{equation}
    \partial_t \phi = \nabla^2 \frac{\delta F}{\delta \phi} + \lambda \nabla^2 \left( \nabla \phi \right)^2 - \gamma \nabla \left[ \nabla^2 \phi \cdot \nabla \phi \right]
\end{equation}
$$

First get rid of powers of gradients with $\left( \nabla \phi \right)^2 \rightarrow \frac{1}{\varepsilon} W(\phi)$.
Then for the last term use product rule and turn $D = \nabla^2 \phi$
into a separate variable. Then we can turn it into $\nabla D \cdot \nabla \phi$. This individual coordinates can then be treated separetly in AMDiS.

$$
\begin{aligned}
    \partial_t \phi &= \nabla^2 K - \gamma \left( D^2 + D^x \partial_x \phi + D^y \phi \right) \\
    K &= M \mu + \frac{\lambda}{\varepsilon} W(\phi) \\
    \mu &= \frac{1}{\varepsilon} W^\prime(\phi) - \varepsilon D \\
    D &= \nabla^2 \phi \\
    D^x &= \partial_x D \\
    D^y &= \partial_y D
\end{aligned}
$$

Here technically not necessary to split $K$ and $\mu$, but it will be once we add the curvature part.

```math
$$
\begin{bmatrix}
    \frac{1}{\tau} + \gamma \left( D^x_n \partial_x + D^y_n \partial_y \right) &
    - \nabla^2 &
    0 & 0 & 0 & 0 \\
    - \frac{\lambda}{\varepsilon} W^\prime (\phi_n) &
    1 &
    - M &
    0 & 0 & 0 \\
    - \frac{1}{\varepsilon} W^{\prime\prime}(\phi_n) &
    0 &
    1 &
    \varepsilon &
    0 & 0 \\
    - \nabla^2 &
    0 & 0 &
    1 &
    0 & 0 \\
    0 & 0 & 0 &
    - \partial_x &
    1 &
    0 \\
    0 & 0 & 0 &
    - \partial_y &
    0 &
    1
\end{bmatrix}
\begin{bmatrix}
    \phi_{n+1} \\
    K_{n+1} \\
    \mu_{n+1} \\
    D_{n+1} \\
    D^x_{n+1} \\
    D^y_{n+1}
\end{bmatrix}
=
\begin{bmatrix}
    \frac{1}{\tau} - \gamma D_n^2 \\
    \frac{\lambda}{\varepsilon} \left[ W(\phi_n) - W^\prime(\phi_n) \phi_n \right] \\
    \frac{1}{\varepsilon} \left[ W^\prime(\phi_n) - W^{\prime\prime}(\phi_n)\phi_n \right] \\
    0 \\
    0 \\
    0
\end{bmatrix}
$$
```

Additionally for the weak form:

$$
\begin{aligned}
    \langle \nabla^2 \phi, v \rangle &= - \langle \nabla \phi, \nabla v \rangle \\
    \langle \partial_x \phi, v \rangle &= - \langle \phi, \partial_x v \rangle
\end{aligned}
$$

for dirichlet and periodic BC.
