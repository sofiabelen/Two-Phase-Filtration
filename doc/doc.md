---
title: Two-dimensional Filtration of Ideal Gas and Pentane
  $\left( C_5 H_{12} \right)$
---

# Physical model

Darcy's law: $$\bm{v}_i = -\frac{1}{\mu_i} \hat K \cdot f_\alpha (s)
    \cdot \nabla P$$ where $\hat K$, the specific permeability. It
depends only on the geometry of the medium. We assume isotropy of space,
so K is a scalar. $\mu$ is the dynamic viscosity.

$i$ - component.

$\alpha$ - phase. (If we had multiple phases, then it would be
$f_\alpha$)

As an approximation, $f_\alpha (s) = s^2$ for the first component, and
$f_\alpha (s) = (1 - s)^2$ for the second.

The continuity equation for each component becomes:
$$\varphi \frac{\partial \rho_i}{\partial t}
    + div (\rho_i \bm{v}_i) = 0$$ where $\rho_i = \frac{m_i}{V}$.

We use the Tait equation to relate liquid density to pressure:
$$\frac{\hat{\rho} - \rho_0}{\hat{\rho}} = C \log_{10}
    \frac{B + P}{B + P_0}$$ where $C = 0.2105$,
$\rho_0 = \frac{1}{67.28 \frac{m^3}{mol}}$, $P_0 = 0.1 MPa$,
$B = 35MPa$, in the case of $C_5H_{12}$.

Ideal gas equation of state:

$$P = \frac{RT}{M} \hat{\rho}$$

# Boundary and Initial Conditions

On the first iteration, we set an initial pressure and molar
composition. Then, we derive the velocities from the pressure gradient
using Darcy's law and the saturation from the densities and molar
composition.

![Boundary Conditions](img/diagram.pdf){#fig:img-diagram-pdf
width="50%"}

## Pressure

**BC:**

$$\begin{cases}
    P = P_{in} &\text{at } y = 0
    \text{ and } x \in [0, 1] \\
    P = P_{out} &\text{at } y = 2 \\
    \frac{\partial P}{\partial x} = 0 &\text{at }
    x = 0, 2 \\
    \frac{\partial P}{\partial y} = 0 &\text{at }
    y = 0 \text{ and } x \in [1, 2]
\end{cases}$$

**IC:**

$$\begin{cases}
    P = P_{out} &\text{at outlet} \\
    P = P_{in}  &\text{at inlet} \\
    P = P_0     &\text{everywhere else}
\end{cases}$$

## Velocities

**IC:** Darcy 2-nd order.

**BC:**

$$\begin{cases}
    u = 0 &\text{at } x = 0, 2 \\
    v = 0 &\text{at } y = 0 \text{ and } x \in [1, 2] \\
    \text{Darcy (2-nd Order FD)}
          &\text{at } y = 2 \\
    \text{Darcy (2-nd Order FD)}
          &\text{at } y = 0
    \text{ and } x \in [0, 1]\\
\end{cases}$$

## Density

We derive the densities from the equations of state.

## Saturation

Boundary condition on the inlet as the molar composition $\psi$:
$$\frac{m_1}{m_2} = \frac{\psi}{1 - \psi}
        \frac{M_1}{M_2}$$ where $M_1$ and $M_2$ represent the molar mass
of each component.

We can derive the densities from the equations of state, then we can
find the saturation.

$$\hat \rho_1 = \frac{m_1}{sV}, \qquad
    \hat \rho_2 = \frac{m_2}{(1 - s)V}$$

$$\frac{\hat \rho_1}{\hat \rho_2} = \frac{m_1}{m_2}
    \frac{1 - s}{s}
    = \frac{\psi M_1}{(1 - \psi) M_2}\frac{1 - s}{s}$$

$$s = \left( 
        \frac{\hat \rho_1 M_2(1 - \psi) }{\hat \rho_2 M_1 \psi}
    + 1 \right)^{-1}$$

**On outlet:**

$$\frac{\partial \alpha}{\partial \bm{n}} = 0$$

Or,

$$\frac{\partial s}{\partial \bm{n}} = 0$$

# Discretization Scheme

We use the finite difference method (FDM).

## Continuity Equation

$$\frac{\partial \rho}{\partial t}
+ \frac{\partial \rho}{\partial x} u
+ \frac{\partial \rho}{\partial y} v
+ \frac{\partial u}{\partial x} \rho
+ \frac{\partial v}{\partial y} \rho = 0$$

$$\varphi \frac{\rho^{n + 1}_{i, j} - \rho^n_{i, j}}{\Delta t}
+ \frac{\rho_{i+1, j}^n - \rho_{i-1,j}^n}{2\Delta x} u_{i,j}
+ \frac{\rho_{i, j+1}^n - \rho_{i,j-1}^n}{2\Delta y} v_{i,j}
+ \frac{u_{i+1, j} - u_{i-1,j}^n}{2\Delta x} \rho_{i,j}^n
+ \frac{v_{i, j+1}^n - v_{i,j-1}^n}{2\Delta y} \rho_{ij}^n = 0$$

## Darcy's Law

$$u_{i,j}^n = -\frac{K}{\mu} f_\alpha (s_{i, j}^n)
    \frac{P_{i + 1, j}^n - P_{i - 1, j}^n}{2\Delta x} \\
$$

$$v_{i,j}^n = -\frac{K}{\mu} f_\alpha(s_{i, j}^n)
    \frac{P_{i, j + 1}^n - P_{i, j - 1}^n}{2\Delta y}$$

### BC

$$v_{i, 1}^n = -\frac{K}{\mu} f_\alpha(s_{i, 1}^n)
    \frac{-3 P_{i, 1}^n + 4 P_{i, 2}^n - P_{i, 3}^n }
    {2 \Delta x} \qquad \text{(Inlet)}$$

$$v_{i, ny}^n = -\frac{K}{\mu} f_\alpha(s_{i, ny}^n)
    \frac{-3 P_{i, ny}^n + 4 P_{i, ny - 1}^n
    - P_{i, ny - 2}^n }{(-2 \Delta x)}
    \qquad \text{(Outlet)}$$

# Finding Pressure using Binary Search

The function *find_pressure* takes as arguments $\rho_1$ and $\rho_2$,
which are defined as

$$\rho_1 = \frac{m_1}{V}, \qquad
\rho_2 = \frac{m_2}{V}
.$$

We try to find the zero of the following function, that takes the
pressure as an argument:

$$f(P) = \frac{\hat{\rho_2} - \rho_0}{\hat{\rho_2}}
    - C \log_{10} \frac{B + P}{B + P_0}
.$$

Here, $\hat{\rho_2} = \frac{m_2}{(1 - s)V}$. In order to determine
$\hat{\rho_2}$, we first find $\hat{\rho_1}$ using the EoS, and from
there, we are able to find the saturation
$s = \frac{\rho_1}{\hat{\rho_1}}$. Lastly, we determine
$\hat{\rho_2} = \frac{\rho_2}{1 - s}$.

# Algorithm

Euler method for discretization with respect to time.

Second order scheme in space, with the use of ghost cells.

1.  Calculate densities for each component using EOS.

2.  Find the pressure and saturation with the help of the Newton Raphson
    method.

3.  Use Darcy's law to calculate velocities.

# Units and Parameters

::: {#tab:label}
                 Temperature                              298 $K$
  ----------------------------------------- -----------------------------------
                  $P_{in}$                               $10^6 Pa$
                  $P_{out}$                              $10^5 Pa$
             Porosity, $\varphi$                            0.7
         Specific Permeability, $K$                     $10^{-12}$
   Dynamic Viscosity of Ideal Gas, $\mu_1$   $1.8 \cdot 10^{-5}$ $Pa \cdot s$
    Dynamic Viscosity of Pentane, $\mu_2$    $2,14 \cdot 10^{-4}$ $Pa \cdot s$
       Molar Mass of Ideal Gas, $M_1$            $0.028$ $\frac{kg}{mol}$
        Molar Mass of Pentane, $M_2$            $0.07215$ $\frac{kg}{mol}$
     Molar Composition at Inlet, $\psi$                    $0.3$

  : Parameters for our simulation.
:::

# TODO

# Conventions

1.  Density:

    $$\rho_i = \frac{m_i}{V}, \qquad
                \hat{\rho_i} = \frac{m_i}{s_i V}
            .$$
