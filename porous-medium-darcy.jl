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
## P = P0 everywhere else
##
## Boundary conditions:
## u = 0 at x = 0, 2 (walls)
## P = Pin at y = 0 and x ∈ [0, 1]
## v = 0 at y = 0 and x ∈ [1, 2]
## P = Pout, dv / dy = 0 at y = 2
## ∂P / ∂y = 0 at y = 0, 2

using DelimitedFiles, PyPlot

function update!(; ρnext, unext, vnext, ρwork, uwork, vwork,
        Δt, Δx, Δy, Pin, Pout, φ, K, μ)
    nx, ny = size(ρnext)
    coeff = 8.314 * 298 / 0.029
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
            ρnext[i, j] = ρc - Δt / φ * (uc * (ρe - ρw) / (2 * Δx) +
                                     vc * (ρn - ρs) / (2 * Δy) +
                                     ρc * (ue - uw) / (2 * Δx) +
                                     ρc * (vn - vs) / (2 * Δy))
        end
    end

    ## Boundary conditions on pressure using ghost cells
    ## P = Pin at y = 0 and x ∈ [0, 1]
    @. @views ρnext[1: nx ÷ 2, 1] = 2 * Pin / coeff - ρnext[1: nx ÷ 2, 2]

    ## P = Pout at y = 2
    @. @views ρnext[:, ny] = 2 * Pout / coeff - ρnext[:, ny - 1]

    ## ∂P / ∂x = 0 at x = 0, 2
    @. @views ρnext[nx, :] = ρnext[nx - 1, :]
    @. @views ρnext[1, :] = ρnext[2, :]

    ## ∂P / ∂y = 0 at y = 0 x ∈ [1, 2]
    @. @views ρnext[nx ÷ 2 + 1: nx, 1] = ρnext[nx ÷ 2 + 1: nx, 2]
    
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
            unext[i, j] = -K / μ * coeff * (ρe - ρw) / (2 * Δx)
    
            vnext[i, j] = -K / μ * coeff * (ρn - ρs) / (2 * Δy)
        end
    end

    ## Boundary conditions on velocity, again using ghost cells
    ## u = 0 at x = 0, 2 (walls)
    @. @views unext[1, :] = -unext[2, :]
    @. @views unext[nx, :] = -unext[nx - 1, :]

    ## v = 0 at y = 0 and x ∈ [1, 2]
    @. @views vnext[1, nx ÷ 2 + 1: nx] = -vnext[2, nx ÷ 2 + 1: nx]

    ## dv / dy = 0 at y = 2
    @views vnext[:, ny] = vnext[:, ny - 1]

    ## dv / dy = 0 at y = 0 and x ∈ [0, 1]
    @views vnext[1: nx ÷ 2, 1] =  vnext[1: nx ÷ 2, 2]
end

function filtration!(; ρ, u, v, Δx, Δy, Δt, nsteps, K, φ, μ, Pin, Pout)
    ## Helper arrays to store values from previous iteration
    ρwork, uwork, vwork = ρ, u, v
    unext = copy(u)
    vnext = copy(v)
    ρnext = copy(ρ)
    
    for t = 1 : nsteps
        update!(; ρnext, unext, vnext, ρwork, uwork, vwork,
                Δt, Δx, Δy, Pin, Pout, φ, K, μ)

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

    nx, ny= size(u)
    
    rgx = range(0, 2, length=nx)
    rgy = range(0, 2, length=ny)
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
    
    savefig("img/porous-medium-darcy.png")
    savefig("img/porous-medium-darcy.svg")
end
let
    ## Time steps
    duration = 2000
    Δt = 0.001
    nsteps = round(Int64, duration / Δt)
    # nsteps = 1000
    
    ## Space grid [0, 2] × [0, 2]
    Δx = 0.05
    Δy = Δx
    nx = convert(Int64, 2 / Δx)
    ny = convert(Int64, 2 / Δy)
    
    ## Porosity
    φ = 0.7
    
    ## Specific permeability
    K = 1e-12
    
    ## Dynamic viscocity (value used: air)
    μ = 18e-6
    
    coeff = 8.314 * 298 / 0.029
    
    ## Pressures at top and bottom, and initial pressure
    Pin = 1e6
    Pout = 1e5
    P0 = Pout
    
    u = zeros(nx, ny)
    v = zeros(nx, ny)
    ρ = zeros(nx, ny)
    
    ρ .= P0 / coeff
    
    
    problem = (
               ρ = ρ,
               u = u,
               v = v,
               Δx = Δx,
               Δy = Δy,
               Δt = Δt,
               nsteps = nsteps,
               K = K,
               μ = μ,
               φ = φ,
               Pin = Pin,
               Pout = Pout
              )
    
    ## Run our solver
    ρ, u, v = filtration!(; problem...)
    
    # writedlm("density.txt", ρ, ' ')

    ## Plotting
    plot(u, v, ρ)
end
