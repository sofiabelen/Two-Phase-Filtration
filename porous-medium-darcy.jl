## Physics:
##
## We simulate 2d flow of ideal gas through a porous medium
## with the help of Darcy's law.
## 
## v⃗ = -μ⁻¹ K̂ ⋅∇P,
##
## where K̂, a second-order tensor in the general case, is
## the specific permeability. It depends only on the
## geometry of the mediumm. Here we assume isotropy
## of space, so K is a scalar. μ is the dynamic viscocity.
##
## Continuity equation:
##
## φ ⋅∂ρ / ∂t + div(ρv⃗) = 0
##
## where φ is the porosity.
##
## So, the equations used are the continuity equation and
## Darcy's equation, which replaces the momentum equation,
## used in our last program.
##
## Another change with respect to the ideal gas program
## is that we will set the boundary conditions using ghost cells.
## This way, we can make use of the central difference scheme for
## boundary conditions which consist of partial derivatives,
## which is second-order.
##
## The boundary conditions consist of walls on the top
## and bottom, and constant pressure on either side.
##
## For now, our velocities will be equal to 0 at the beginning.
##
## ---------
## ->     ->
## ->     ->
## ---------
##
## Initial conditions:
## v⃗ = (u, v) = 0⃗ everywhere
## P = Pin, Pout at x = 0, 2 (respectively)
## P = P0 everywhere else
##
## Boundary conditions:
## ∂u / ∂x = 0 at x = 0, 2
## P = Pin at x = 0
## P = Pout at x = 2
## u = 0 at y = 0, 2
## ∂P / ∂y = 0 at y = 0, 2
## v = 0 at x = 0, 2 and y = 0, 2 (is this okay???)

using PyPlot

## Time steps
duration = 10.0
Δt = 0.001
nsteps = convert(Int64, duration / Δt)

## Space grid [0, 2] × [0, 2]
Δx = 0.05
Δy = Δx
nx = convert(Int64, 2.0 / Δx)
ny = convert(Int64, 2.0 / Δy)

## Porosity
φ = 0.7

## Specific permeability
K = 1e-17

## Dynamic viscocity (value used: air)
μ = 0.018

## Pressures at left and right side.
Pin = 1.5
Pout = 0.5
P0 = 1.0

u = zeros(nx, ny)
v = zeros(nx, ny)
ρ = zeros(nx, ny)

for i = 2 : nx - 1
    ρ[i, :] .= Pin + (Pout - Pin) * i / nx
end

function update!(ρnext, unext, vnext, ρwork, uwork, vwork,
        Δt, Δx, Δy, nx, ny)
    for j = 2 : ny - 1
        for i = 2 : nx - 1
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
            ρw = ρwork[i - 1, j]
            ρn = ρwork[i, j + 1]
            ρs = ρwork[i, j - 1]
    
            uc = uwork[i, j]
            ue = uwork[i + 1, j]
            uw = uwork[i - 1, j]
    
            vc = vwork[i, j]
            vn = vwork[i, j + 1]
            vs = vwork[i, j - 1]
    
            ## Continuity equation to calculate ρ
            ## φ ⋅∂ρ / ∂t + div(ρv⃗) = 0
            ρnext[i, j] = ρc - Δt / φ * (uc * (ρe - ρw) / (2.0 * Δx) +
                                     vc * (ρn - ρs) / (2.0 * Δy) +
                                     ρc * (ue - uw) / (2.0 * Δx) +
                                     ρc * (vn - vs) / (2.0 * Δy))
        end
    end
    
    ## Question about view: https://stackoverflow.com/questions/65979923/what-is-happening-behind-the-scenes-in-julias-view-function-a3-viewa

    ## Boundary conditions using ghost cells
    ## ρ = ρleft, ρright at x = 0, 2
    ρnext[1, :] .= 2.0 * Pin .- ρnext[2, :]
    ρnext[nx, :] .= 2.0 * Pout .- ρ[nx - 1, :]
    ## ∂ρ / ∂y = 0 at y = 0, 2
    ρnext[:, ny] .= ρnext[:, ny - 2]
    ρnext[:, 1] .= ρnext[:, 3]
    
    ## Here we use a central difference scheme to discretize the second
    ## order derivatives, so we iterate over the 'inner' grid
    for j = 2 : ny - 1
        for i = 2 : nx - 1
            ## Auxiliar variables
            ρe = ρwork[i + 1, j]
            ρw = ρwork[i - 1, j]
            ρn = ρwork[i, j + 1]
            ρs = ρwork[i, j - 1]
    
            ## Darcy's law
            ## v⃗ = -μ⁻¹ K̂ ⋅∇P
            unext[i, j] = -K / μ * (ρe - ρw) / (2 * Δx)
    
            vnext[i, j] = -K / μ * (ρn - ρs) / (2 * Δy)
        end
    end

    ## Boundary conditions on u, again using ghost cells
    ## ∂u / ∂x = 0 at x = 0, 2
    unext[1, :] = unext[3, :]
    unext[nx, :] = unext[nx - 2, :]
    ## u = 0 at y = 0, 2
    u[:, ny] = -u[:, ny - 1]
    u[:, 1] = -u[:, 2]
end

function filtration!(ρ, u, v, Δx, Δy, Δt, nsteps, nx, ny)
    ## Helper arrays to store values from previous iteration
    ρwork, uwork, vwork = ρ, u, v
    unext = copy(u)
    vnext = copy(v)
    ρnext = copy(ρ)
    
    for t = 1 : nsteps
        update!(ρnext, unext, vnext, ρwork, uwork, vwork,
                Δt, Δx, Δy, nx, ny)

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
    quiver(x, y, u', v')
    
    ## Streamplot of velocity vector field
    #streamplot(x, y, u', v')
    
    ## Density contour map
    pos = ax.contourf(X, Y, ρ', alpha=0.5, cmap=matplotlib.cm.viridis)
    fig.colorbar(pos, ax=ax)
    cp = contour(X, Y, ρ', cmap=matplotlib.cm.viridis)

    PyPlot.title("2d filtration of ideal gas through a porous medium")

    # PyPlot.svg(true)
    # savefig("plot.svg")
    
    savefig("img/porous-medium-darcy.png")
    savefig("img/porous-medium-darcy.svg")
end

## Run our solver
ρ, u, v = filtration!(ρ, u, v, Δx, Δy, Δt, nsteps, nx, ny)

## Plotting
plot(u, v, ρ)
