### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 04bbf20d-ebb9-4acc-b821-27ee3d5b3b0a
md"""
# Boundary Conditions Documentation

We want to create a general way to update the BC for readability and in order to minimize errors.

``
\left. \left( af(x, y) + b\frac{\partial f}{\partial n} \right) \right\vert_{(x, y) \in \partial \Omega} = g
``

#### Discretization

``
a \frac{f(\vec{r} - \vec{h})  + f(\vec{r} + \vec{h})}{2} +
b \frac{f(\vec{r} + \vec{h}) - f(\vec{r} - \vec{h})}{2h} = g
``

##### Example: Right Boundary

``
\vec{h} = (\frac{\Delta x}{2}, 0)
``

``
a \frac{f_{nx} + f_{nx - 1}}{2} + b \frac{f_{nx} - f_{nx - 1}}{\Delta x} = g
``

``
f_{nx} = \frac{2 \Delta x g + f_{nx - 1} (2b - a \Delta x)}{a \Delta x + 2b}
``

##### General Formula

``
f_{ghost} = \frac{2 \Delta h g + f_{inside} (2b - a \Delta h)}{a \Delta h + 2b}
``

#### Implementation

The struct Boundary will hold 4 arrays, one for each border: left, right, top, bottom, in the form of a 2D array.

x[1] -> left

x[2] -> right

x[3] -> top

x[4] -> bottom

Then, the struct BoundaryCondition will hold the parameters a, b and g, each of type Boundary.
"""

# ╔═╡ 1f5ca81e-3d74-11ec-2394-2d82e6fcccbd
struct Boundary{T}
    x::Array{T, 2}
end

# ╔═╡ 50013c2c-befe-4483-bc6f-764a31a33ccd
struct BoundaryCondition{T}
    a::Boundary{T}
    b::Boundary{T}
end

# ╔═╡ 923b99aa-52dc-4068-b888-b5cc9d8c516a
md"""
#### Algorithm

1. Enforce the BC: ``\frac{\partial P}{\partial n} = 0`` on the walls.

2. Obtain ``\hat\rho_{\alpha, \partial \Omega} = \hat\rho_{\alpha, \partial \Omega}(T, P_{\partial \Omega}) `` from their respective equations of state. ``\alpha = 1, 2`` - component. On the inlet and outlet, the initial pressure is given as ``P_{in}`` and ``P_{out}``, respectively. On the walls, we derived it in the previous step.

3. Obtain the boundary saturation.

On the inlet:

``
s_{in} = s_{in}(\psi_{in},
\hat \rho_{1, in}, \hat \rho_{2, in})
= \left( 
        \frac{\hat \rho_{1, in} M_2(1 - \psi_{in}) }{\hat \rho_{2, in} M_1 \psi_{in}}
    + 1 \right)^{-1}
``

``
\left. \left( 1 \cdot s +
0 \cdot \frac{\partial s}{\partial n} \right) \right\vert_{inlet}
= s_{in} 
``

On the outlet and walls:

``
\frac{\partial s}{\partial n} = 0
``

``
\left. \left( 0 \cdot s +
1 \cdot \frac{\partial s}{\partial n} \right) \right\vert_{outlet, walls}
= 0
``

4. Find ``\rho_{\alpha, \partial \Omega} = \hat\rho_{\alpha, \partial \Omega} s_\alpha``.

``
\left. \left( 1 \cdot \rho_\alpha +
0 \cdot \frac{\partial \rho_\alpha}{\partial n} \right) \right\vert_{inlet}
= \rho_{\alpha, \partial \Omega}
``

5. Update densities at boundary.

``
\rho_{k, \partial \Omega} = \hat \rho_{k, \partial \Omega} \cdot
s_{k, \partial \Omega}
``

6. Find ``P_{\partial \Omega}`` on the inlet and outlet from the densities at the boundary (using Newton-Raphson).

7. Enforce BC for velocity using Darcy's Law on inlet and outlet

``
    v_{in} = -\frac{K}{\mu} f_\alpha (s_{in})
    \frac{P_{i, 2} - P_{i, 1}^n}{\Delta y}
``

``
    v_{out} = -\frac{K}{\mu} f_\alpha (s_{out})
    \frac{P_{i, ny} - P_{i, ny - 1}^n}{\Delta y}
``

``
\left. \left( 1 \cdot v +
0 \cdot \frac{\partial v}{\partial n} \right) \right\vert_{in, out}
= v_{in, out}
``

"""

# ╔═╡ 2ef6d608-3aef-4647-a523-47fab6ff1154
md"""
## Todo
1. Now it runs... but nothing changes with time :c
"""

# ╔═╡ f732b89b-5361-4dc5-9051-62fa8ae698cb


# ╔═╡ c427e959-cc5e-4eeb-b50c-8c06aec042e8


# ╔═╡ 6f6caa2a-f2ed-4fa0-9964-0af022384749


# ╔═╡ ce3afddf-9b25-4468-8777-b77a6898da5f


# ╔═╡ Cell order:
# ╟─04bbf20d-ebb9-4acc-b821-27ee3d5b3b0a
# ╠═1f5ca81e-3d74-11ec-2394-2d82e6fcccbd
# ╠═50013c2c-befe-4483-bc6f-764a31a33ccd
# ╟─923b99aa-52dc-4068-b888-b5cc9d8c516a
# ╟─2ef6d608-3aef-4647-a523-47fab6ff1154
# ╠═f732b89b-5361-4dc5-9051-62fa8ae698cb
# ╠═c427e959-cc5e-4eeb-b50c-8c06aec042e8
# ╠═6f6caa2a-f2ed-4fa0-9964-0af022384749
# ╠═ce3afddf-9b25-4468-8777-b77a6898da5f
