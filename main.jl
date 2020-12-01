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
duration = 5.0
Δt = 0.001
nsteps = convert(Int64, duration / Δt)

## Space grid [0, 2] × [0, 2]
Δx = 0.05
Δy = Δx
nx = convert(Int64, 2.0 / Δx)
ny = convert(Int64, 2.0 / Δy)

## Viscocity
ν = 0.1

## Atmospheric pressure
p0 = 1.0

## Velocity at the top
u0 = 0.05

## Initial conditions:
## u = u0 = 0.05 at y = 2
## v⃗ = (u, v) = 0⃗ everywhere else
## ρ = p = p0 everywhere (Atmospheric pressure?)
##
## Boundary conditions:
## u = 1 at y = 2
## ∂p / ∂y = 0 at y = 0
## ρ = p = p0 at y = 2
## ∂p / ∂x = 0 at x = 0, 2

u = zeros(nx, ny)
v = zeros(nx, ny)
ρ = zeros(nx, ny)

for i = 1 : nx
    for j = 1: ny
        ρ[i, j] = p0
    end
end

for i = 1 : nx
    u[i, ny] = u0
end

function update!(ρnext, unext, vnext, ρwork, uwork, vwork,
        Δt, Δx, Δy, u0, ν, nx, ny)
    for i = 1 : nx - 1
        for j = 1 : ny - 1
            ## Auxiliary variables to avoid typos when
            ## working with partial derivatives
            ##
            ## NW--N--NE
            ##     |
            ## W---C---E
            ##     |
            ## SW--S--SE
    
            ρc = ρwork[i, j]
            ρe = ρwork[i + 1, j]
            ρn = ρwork[i, j + 1]
    
            uc = uwork[i, j]
            ue = uwork[i + 1, j]
            un = uwork[i, j + 1]
    
            vc = vwork[i, j]
            ve = vwork[i + 1, j]
            vn = vwork[i, j + 1]
    
            ## Continuity equation to calculate ρ
            ## ∂ρ / ∂t + div(ρv⃗) = 0
            ρnext[i, j] = ρc - Δt * (uc * (ρe - ρc) / Δx +
                                     vc * (ρn - ρc) / Δy +
                                     ρc * (ue - uc) / Δx +
                                     ρc * (vn - vc) / Δy)

        end
    end
    
    ## Enforcing boundary conditions for ρ
    for j = 1 : ny
        ρnext[1, j] = ρnext[2, j]
        ρnext[nx, j] = ρnext[nx - 1, j]
    end
    
    for i = 1 : nx
        ρnext[i, 1] = ρnext[i, 2]
        ρnext[i, ny] = p0
    end
    
    ## Here we use a central difference scheme to discretize the second
    ## order derivatives, so we iterate over the 'inner' grid
    for i = 2 : nx - 1
        for j = 2 : ny - 1
            ## Auxiliar variables
            ρc = ρwork[i, j]
            ρe = ρwork[i + 1, j]
            ρn = ρwork[i, j + 1]
    
            uc = uwork[i, j]
            ue = uwork[i + 1, j]
            uw = uwork[i - 1, j]
            un = uwork[i, j + 1]
            us = uwork[i, j - 1]
    
            vc = vwork[i, j]
            ve = vwork[i + 1, j]
            vw = vwork[i - 1, j]
            vn = vwork[i, j + 1]
            vs = vwork[i, j - 1]
    
            ## Momentum equations 
            ## ∂v⃗ / ∂t + (v⃗ ⋅∇)v⃗ = - (1 / ρ) ∇P + ν∇²v⃗
            viscous_term_u = ν * ((uw - 2.0 * uc + ue) / Δx^2 +
                                  (un - 2.0 * uc + us) / Δy^2)
    
            viscous_term_v = ν * ((vw - 2.0 * vc + ve) / Δx^2 +
                                  (vn - 2.0 * vc + vs) / Δy^2)
    
            pressure_term_u = -(1 / ρc) * (ρe - ρc) / Δx
    
            pressure_term_v = -(1 / ρc) * (ρn - ρc) / Δy
    
            inertia_term_u = uc * (ue - uc) / Δx +
                             vc * (un - uc) / Δy
    
            inertia_term_v = uc * (ve - vc) / Δx +
                             vc * (vn - vc) / Δy
    
            unext[i, j] = uc + Δt * (-inertia_term_u + pressure_term_u +
                                     viscous_term_u)
    
            vnext[i, j] = vc + Δt * (-inertia_term_v + pressure_term_v +
                                     viscous_term_v)
        end
    end
end

function cavity_flow!(ρ, u, v, Δx, Δy, Δt, nsteps, nx, ny, u0, ν)
    ## Helper arrays to store values from previous iteration
    ρwork, uwork, vwork = ρ, u, v
    unext = copy(u)
    vnext = copy(v)
    ρnext = copy(ρ)
    
    for t = 1 : nsteps
        update!(ρnext, unext, vnext, ρwork, uwork, vwork,
                Δt, Δx, Δy, u0, ν, nx, ny)

        ## Swap next and current and continue iteration
        ρnext, ρwork = ρwork, ρnext
        unext, uwork = uwork, unext
        vnext, vwork = vwork, vnext
    end

    ## If in the end ρwork is not the original ρ,
    ## copy the most recent data into ρ, u, v
    if ρwork !== ρ
        ρ .= ρwork
        u .= uwork
        v .= vwork
    end

    return ρ, u, v
end

function plot(u, v, ρ)
    fig = PyPlot.figure(figsize=(10, 10))
    ax = PyPlot.axes()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    #legend(loc="upper right",fancybox="true")
    
    rgx = range(0, 2.0, length=nx)
    rgy = range(0, 2.0, length=ny)
    x = [rgx;]'
    y = [rgy;]'
    X = repeat([rgx;]', length(x))
    Y = repeat([rgy;],1, length(y))
    
    ## Velocity vector field
    ## scale = 100.0
    ## quiver(x, y, u', v')
    
    ## Streamplot of velocity vector field
    streamplot(x, y, u', v')
    contour(X, Y, ρ')
    
    ## Density contour map
    pos = ax.contourf(X, Y, ρ', alpha=0.5, cmap=matplotlib.cm.viridis)
    fig.colorbar(pos, ax=ax)
    cp = contour(X, Y, ρ', cmap=matplotlib.cm.viridis)

    PyPlot.title("Cavity Flow for Ideal Gas")

    # PyPlot.svg(true)
    # savefig("plot.svg")
    
    savefig("plot.png")
end

## Run our solver
ρ, u, v = cavity_flow!(ρ, u, v, Δx, Δy, Δt, nsteps, nx, ny, u0, ν)

## Plotting
plot(u, v, ρ)
