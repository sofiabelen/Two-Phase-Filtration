## Physics:
##
## We simulate 2d cavity flow for compressible fluid,
## more precisely, for an ideal gas at constant temperature.
## Therefore, from the equation of state, we know that
## P / ρ = 1 (if we choose the units such that the constants
## terms are reduced to 1). 
## 
## The equations used are the continuity and the momentum equations
## First order discretization schemes are used.
##
## The velocity in the x-direction on the top is constant,
## there are walls on the bottom, left and right edges.
##
## →→→→→→→→
## |      |
## |      |
## --------

using PyPlot

## Time steps
duration = 1.0
Δt = 0.001
nt = convert(Int64, duration / Δt)

## Space grid [0, 2] × [0, 2]
Δx = 0.1
Δy = Δx
nx = convert(Int64, 2.0 / Δx)
ny = convert(Int64, 2.0 / Δy)

## Viscocity
ν = 0.1

## Atmospheric pressure
p0 = 1.0

## Initial conditions:
## u = 1 at y = 2
## v⃗ = (u, v) = 0⃗ everywhere else
## ρ = p = p0 everywhere (Atmospheric pressure?)
u = zeros(nx, ny)
v = zeros(nx, ny)
ρ = zeros(nx, ny)

for i = 1 : nx
    for j = 1: ny
        ρ[i, j] = p0
    end
end

for i = 1 : nx
    u[i, ny] = 1.0
end

## Boundary conditions:
## u = 1 at y = 2
## ∂p / ∂y = 0 at y = 0
## ρ = p = p0 at y = 2
## ∂p / ∂x = 0 at x = 0, 2

## Helper arrays to store values from previous iteration
un = u
vn = v
ρn = ρ

## Solver
for t = 1 : nt
    for i = 1 : nx - 1
        for j = 1 : ny - 1
            ## Continuity equation to calculate ρ
            ## ∂ρ / dt + div(ρv⃗) = 0
            global ρ[i, j] = ρn[i, j] - Δt * ((ρn[i + 1, j] - ρn[i, j]) / Δx *
                        un[i, j] + (ρn[i, j + 1] - ρn[i, j]) / Δy +
                        (un[i + 1, j] - un[i, j]) / Δx * ρn[i, j] +
                        (vn[i, j + 1] - vn[j, j]) / Δy * ρn[i, j])
        end
    end

    ## Enforcing boundary conditions for ρ
    # for j = 1 : ny
    #     global ρ[1, j] = ρ[2, j]
    #     global ρ[nx, j] = ρ[nx - 1, j]
    # end

    for i = 1 : nx
        global ρ[i, 1] = ρ[i, 2]
        global ρ[i, ny] = p0
    end
    
    ## Here we use a central difference scheme to discretize the second
    ## order derivatives, so we iterate over the 'inner' grid
    for i = 2 : nx - 1
        for j = 2 : ny - 1
            ## Momentum equations 
            ## ∂v⃗ / ∂t + (v⃗ ⋅∇)v⃗ = - (1 / ρ) ∇P + ν∇²v⃗
            global u[i, j] = un[i, j] + Δt * (-un[i, j] * (un[i + 1, j] - un[i, j]) / Δx -
                        vn[i, j] * (un[i, j + 1] - un[i,j]) / Δy -
                        1 / ρn[i, j] * (ρn[i + 1, j] - ρn[i, j]) / Δx +
                        ν * ((un[i + 1, j] - 2.0 * un[i, j] + un[i - 1, j]) / Δx^2 +
                             (un[i, j + 1] - 2.0 * un[i, j] + un[i, j - 1]) / Δy^2))

            global v[i, j] = vn[j, j] + Δt * (-un[i, j] * (vn[i + 1, j] - vn[i, j]) / Δy -
                        vn[i, j] * (vn[i, j + 1] - vn[i, j]) / Δy -
                        1 / ρn[i, j] * (ρn[i, j + 1] - ρn[i, j]) / Δy +
                        ν * ((vn[i + 1, j] - 2.0 * vn[i, j] + vn[i - 1, j]) / Δx
                             + (vn[i, j + 1] - 2.0 * vn[i, j] + vn[i, j - 1]) / Δy))
        end
    end

    ## Enforcing boundary conditions for v⃗
    for i = 1 : nx
        u[i, ny] = 1.0
    end

    global un = u
    global vn = v
    global ρn = ρ
end

## Plotting
fig = figure(figsize=(10, 10))
ax = PyPlot.axes()
ax.set_xlabel("x")
ax.set_ylabel("y")
legend(loc="upper right",fancybox="true")

rgx = range(0, 2.0, length=nx)
rgy = range(0, 2.0, length=ny)
x = [rgx;]'
y = [rgy;]'
X = repeat([rgx;]', length(x))
Y = repeat([rgy;],1, length(y))

## Velocity vector field
# scale = 100.0
quiver(x, y, u', v')

## Density contour map
cp = contour(X, Y, ρ')

# PyPlot.svg(true)
# savefig("plot.svg")

savefig("plot.png")
