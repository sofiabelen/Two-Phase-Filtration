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

mutable struct System{T<:Float64}
    u::Array{T, 3}
    v::Array{T, 3}
    ρ::Array{T, 3}
    P::Array{T, 2}
    s::Array{T, 2}
end

Base.copy(sys::System) = System(sys.u, sys.v, sys.ρ, sys.P, sys.s)

function binarySearch(; f, left, right, niter=50, eps=1e-6)
    mid = left
    for i = 1 : niter
        mid = left + (right - left) / 2

        if abs((right - left) / mid) < eps
            return mid
        elseif f(mid) < 0
            left = mid
        elseif f(mid) > 0
            right = mid
        else
            return mid
        end
    end
    return mid
end

# function findPressure(; m₁::Float64, m₂::Float64, V::Float64,
#         left=500, right=1000, niter=50, eps=1e-6)
#     ## We want to find "any" zero of this function
#     function f(ρ₂::Float64)
#         ρ₁ = m₁ * ρ₂ / (V * ρ₂ - m₂)
#         P = 8.314 * 298 / 0.029 * ρ₁
#         ρ₀ = 616.18
#         P₀ = 1e5
#         return (ρ₂ - ρ₀) / ρ₂ - 0.2105 * log10((35e6 + P) / (35e6 + P₀))
#     end
# 
#     ρ₂ = binarySearch(; f, left, right, eps)
#     ρ₁ = m₁ * ρ₂ / (V * x₀ - m₂)
#     P = 8.314 * 298 / 0.029 * ρ₁
#     return P
# end

function findPressure(ρ̂₁::Float64, ρ̂₂::Float64, V::Float64,
        left=1e4, right=1e7, niter=50, eps=1e-6)
    function f(P::Float64)
        ρ₁ = 0.029 / 8.314 / 298 * P
        s = ρ̂₁ / ρ₁
        ρ₂ = ρ̂₂ / (1 - s)
        ρ₀ = 616.18
        P₀ = 1e5
        return (ρ₂ - ρ₀) / ρ₂ - 0.2105 * log10((35e6 + P) / (35e6 + P₀))
    end

    P = binarySearch(; f, left, right, eps)
    ρ₁ = 0.029 / 8.314 / 298 * P
    s = ρ̂₁ / ρ₁
    return P, s
end

function findDensityofLiquid(P::Float64, s::Float64,
        left=0, right=1000, niter=50, eps=1.e-6)
    function f(ρ₂)
        ρ₀ = 616.18
        P₀ = 1e5
        return (ρ₂ - ρ₀) / ρ₂ - 0.2105
            * log10((35e6 + P) / (35e6 + P₀))
    end

    ρ₂ = binarySearch(; f, left, right, eps)
    ρ̂₂ = (1 - s) * ρ₂
    return ρ̂₂
end

function cnsweRename(arr::Array{Float64, 2}, i::Int, j::Int)
    ## Avoid typos when working with partial derivatives
    ##
    ## NW--N--NE
    ##     |
    ## W---C---E
    ##     |
    ## SW--S--SE
    
    c = arr[i, j]
    n = arr[i, j + 1]
    s = arr[i, j - 1]
    w = arr[i - 1, j]
    e = arr[i + 1, j]
    return c, n, s, w, e
end

#### Question: here ρ represent ρ = m / V? or (vs ρ = m / (sV))?
## Continuity equation to calculate ρ
## φ ⋅∂ρ / ∂t + div(ρv⃗) = 0
function continuityEquation!(Δt, Δx, Δy, φ, ρnext::Array{Float64, 2},
        ρwork::Array{Float64, 2}, uwork::Array{Float64, 2},
        vwork::Array{Float64, 2})
    nx, ny = size(ρnext)
    for j = 2 : ny - 1
        for i = 2 : nx - 1
            ρc, ρn, ρs, ρw, ρe = cnsweRename(ρwork, i, j)
            uc, un, us, uw, ue = cnsweRename(uwork, i, j)
            vc, vn, vs, vw, ve = cnsweRename(vwork, i, j)
    
            ρnext[i, j] = ρc - Δt / φ * (uc * (ρe - ρw) / (2 * Δx) +
                                     vc * (ρn - ρs) / (2 * Δy) +
                                     ρc * (ue - uw) / (2 * Δx) +
                                     ρc * (vn - vs) / (2 * Δy))
        end
    end
end

## Darcy's law
## v⃗ = -μ⁻¹ K̂ ⋅∇P
function darcy(f, unext::Array{Float64, 2}, vnext::Array{Float64, 2},
        Pnext::Array{Float64, 2}, s::Array{Float64, 2},
        K::Float64, μ::Float64, Δx::Float64, Δy::Float64)
    nx, ny = size(unext)
    for j = 2 : ny - 1
        for i = 2 : nx - 1
            ## можно через (ρwork + ρnext) / 2
            Pc, Pn, Ps, Pw, Pe = cnsweRename(Pnext, i, j)

            unext[i, j] = -K / μ * f(s[i, j])
                * (Pe - Pw) / (2 * Δx)
            vnext[i, j] = -K / μ * f(s[i, j])
                * (Pn - Ps) / (2 * Δy)
        end
    end
end

function boundaryDensity!(ρnext::Array{Float64, 2})
    nx, ny = size(ρnext)
    ## ρ = 0 'behind walls'
    @. @views ρnext[1, :] = 0
    @. @views ρnext[nx, :] = 0
    @. @views ρnext[nx ÷ 2 + 1: nx, 1] = 0

    ## ∂ρ / ∂y = 0 at entrance and exit
    @. @views ρnext[:, ny] = ρnext[:, ny - 1]
    @. @views ρnext[1: nx ÷ 2, 1] = ρnext[1: nx ÷ 2, 2]
end

function boundaryPressure(Pnext::Array{Float64, 2},
        Pin::Float64, Pout::Float64)
    nx, ny = size(Pnext)

    ## P = Pin at y = 0 and x ∈ [0, 1]
    @. @views Pnext[1: nx ÷ 2, 1] =
        2 * Pin - Pnext[1: nx ÷ 2, 2]

    ## P = Pout at y = 2
    @. @views Pnext[:, ny] =
        2 * Pout - Pnext[:, ny - 1]

    ## ∂P / ∂x = 0 at x = 0, 2
    @. @views Pnext[nx, :] = Pnext[nx - 1, :]
    @. @views Pnext[1, :] = Pnext[2, :]

    ## ∂P / ∂y = 0 at y = 0 x ∈ [1, 2]
    @. @views Pnext[nx ÷ 2 + 1: nx, 1] =
        Pnext[nx ÷ 2 + 1: nx, 2]
end

function boundaryVelocity(unext::Array{Float64, 2},
        vnext::Array{Float64, 2})
    nx, ny = size(unext)
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

function update!(ρnext, unext, vnext, Pnext, snext,
        ρwork, uwork, vwork, Pwork, swork,
        Δt, Δx, Δy, Pin, Pout, φ, K, μ)
    nx, ny = size(ρnext)
    coeff = 8.314 * 298 / 0.029

    #### Question: How to vectorize and get rid of for loop??
    for k = 1 : 2
        continuityEquation!(Δt, Δx, Δy, φ,
                            ρnext[k, :, :], ρwork[k, :, :],
                            uwork[k, :, :], vwork[k, :, :])
    end
    for k = 1 : 2
        boundaryDensity!(ρnext[k, :, :])
    end
    ## Question: how to vectorize this?
    ## @. Pnext, snext = findPressure(ρnext[1, :, :], ρnext[2, :, :], Δx * Δy)
    for j = 1 : ny
        for i = 1 : nx
            Pnext[i, j], snext[i, j] =
                findPressure(ρnext[1, i, j],
                             ρnext[2, i, j], Δx * Δy)
        end
    end
    boundaryPressure(Pnext, Pin, Pout)
    
    f1(s) = s^2
    f2(s) = (1 - s)^2
    f = [f1, f2]
    for k = 1 : 2
        darcy(f[k], unext[k, :, :], vnext[k, :, :], Pnext,
              snext, K, μ, Δx, Δy)
    end
    for k = 1 : 2
        boundaryVelocity(unext[k, :, :], vnext[k, :, :])
    end
end

function filtration!(; Δx, Δy, Δt, nsteps, K, φ, μ,
        Pin, Pout, sys)
    ## Helper arrays to store values from previous iteration
    syswork = sys
    sysnext = copy(sys)
    
    for t = 1 : nsteps
        update!(sysnext.ρ, sysnext.u, sysnext.v, sysnext.P,
                sysnext.s, syswork.ρ, syswork.u, syswork.v,
                syswork.P, syswork.s,
                Δt, Δx, Δy, Pin, Pout, φ, K, μ)

        ## Swap next and work and continue iteration
        sysnext, syswork = syswork, sysnext
    end

    ## If in the end syswork is not the original sys,
    ## copy the most recent data into sys
    if syswork !== sys
        sys = syswork
    end
end

function plot(u::Array{Float64, 3}, v::Array{Float64, 3},
        P::Array{Float64, 2})
    nx, ny = size(u[1, :, :])

    fig = PyPlot.figure(figsize=(10, 10))
    ax = PyPlot.axes()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    
    rgx = range(0, 2, length=nx)
    rgy = range(0, 2, length=ny)
    x = [rgx;]'
    y = [rgy;]'
    X = repeat([rgx;]', length(x))
    Y = repeat([rgy;],1, length(y))
    
    ## Velocity vector field
    ## scale = 100.0
    quiver(x, y, u[1, :, :]', v[1, :, :]', color="r")
    quiver(x, y, u[2, :, :]', v[2, :, :]', color="b")
    
    ## Pressure contour map
    pos = ax.contourf(X, Y, P', alpha=0.5,
        cmap=matplotlib.cm.viridis)
    fig.colorbar(pos, ax=ax)
    cp = contour(X, Y, P', cmap=matplotlib.cm.viridis)

    PyPlot.title("2 Phase Filtration: Liquid and Ideal Gas")
    
    savefig("img/2phase-filtration.png", dpi=200)
    savefig("img/2phase-filtration.svg")
end

let
    duration = 2000
    Δt = 0.001
    nsteps = round(Int64, duration / Δt)
    nsteps = 1
    
    ## Space grid [0, 2] × [0, 2]
    Δx = 0.05
    Δy = Δx
    nx = round(Int64, 2 / Δx)
    ny = round(Int64, 2 / Δy)
    
    ## Porosity
    φ = 0.7
    
    ## Specific permeability
    K = 1e-12
    
    ## Dynamic viscocity (for air)
    μ = 18e-6
    
    ## Pressures at top and bottom, and initial pressure
    Pin = 1e6
    Pout = 1e5
    P0 = (Pin - Pout) / 2

    ## Initial saturation
    s0 = 0.75

    coeff = 8.314 * 298 / 0.029

    u = zeros(2, nx, ny)
    v = copy(u)
    P = fill(P0, nx, ny)
    s = fill(s0, nx, ny)
    ρ = fill(P0 / coeff, 2, nx, ny)
    @. ρ[2, :, :] = findDensityofLiquid(P, s)
    boundaryPressure(P, Pin, Pout)

    sys = System(u, v, ρ, P, s)

    Parameters = (
                  Δt = Δt,
                  nsteps = nsteps,
                  Δx = Δx,
                  Δy = Δy,
                  φ = φ,
                  K = K,
                  μ = μ,
                  Pin = Pin,
                  Pout = Pout,
                  sys = sys
                 )

    filtration!(; Parameters...)
    writedlm("pressure.txt", sys.P, ' ')
    writedlm("u1.txt", sys.u[1, :, :], ' ')
    writedlm("u2.txt", sys.u[2, :, :], ' ')
    writedlm("density1.txt", sys.ρ[1, :, :], ' ')
    writedlm("density2.txt", sys.ρ[2, :, :], ' ')
    writedlm("s.txt", sys.s, ' ')
    plot(sys.u, sys.v, sys.s)
end
